#include "cldera_profiling_interface.hpp"

#include "cldera_profiling_context.hpp"
#include "cldera_profiling_session.hpp"
#include "cldera_profiling_archive.hpp"
#include "cldera_pathway_factory.hpp"
#include "stats/cldera_register_stats.hpp"

#include "timing/cldera_timing_session.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_session.hpp>
#include <ekat/ekat_assert.hpp>

#include <cstring>
#include <fstream>

extern "C" {

namespace cldera {

inline ProfilingSession& get_session () {
  return ProfilingSession::instance();
}

void cldera_switch_context_c (const char*& name) {
  get_session().switch_context(name);
}

ProfilingContext& get_curr_context () {
  return get_session().get_curr_context();
}

void cldera_init_c (const char*& context_name,
                    const MPI_Fint fcomm,
                    const int case_t0_ymd, const int case_t0_tod,
                    const int run_t0_ymd, const int run_t0_tod,
                    const int stop_ymd, const int stop_tod)
{
  // Convert F90 comm to C comm and create ekat::Comm
  MPI_Comm mpiComm = MPI_Comm_f2c(fcomm);
  EKAT_REQUIRE_MSG (mpiComm!=nullptr, "Error! Input fortran comm seems invalid.\n");
  ekat::Comm comm(mpiComm);

  if (comm.am_i_root()) {
    printf(" [CLDERA] Initializing profiling context '%s' ...\n",context_name);
  }

  // Init session (if not already done)
  auto& s = get_session();
  s.init_session(comm);

  // Add context right away, even before checking if an input file is present.
  // This way, we can always check if this context is inited or not (we can't
  // do that if we don't create a context)
  auto& c = s.add_context(context_name);

  //TODO: make the filename configurable
  std::string filename = "./cldera_profiling_config.yaml";
  if (not std::ifstream(filename).good()) {
    if (comm.am_i_root()) {
      printf(" [CLDERA] WARNING: could not open './cldera_profiling_config.yaml'.\n"
             "   -> Profiling will do nothing.\n");
    }
    return;
  }
  auto session_params = ekat::parse_yaml_file(filename);
  const auto& context_params_filename = session_params.get<std::string>(context_name);
  if (not std::ifstream(context_params_filename).good()) {
    if (comm.am_i_root()) {
      printf(" [CLDERA] WARNING: could not open '%s'.\n"
             "   -> Profiling will do nothing for context '%s'.\n",
             context_params_filename.c_str(),context_name);
    }
    return;
  }
  auto params = ekat::parse_yaml_file(context_params_filename);

  c.init(comm,params);
  std::string timer_name = context_name;
  timer_name += "::init";
  c.timing().start_timer(timer_name);

  using stat_ptr_t = std::shared_ptr<FieldStat>;
  using requests_t = std::map<std::string,std::vector<stat_ptr_t>>;
  using vos_t = std::vector<std::string>;

  TimeStamp case_t0 (case_t0_ymd,case_t0_tod);
  TimeStamp run_t0 (run_t0_ymd,run_t0_tod);
  TimeStamp stop (stop_ymd,stop_tod);
  c.create<TimeStamp>("run_t0",run_t0);
  c.create<TimeStamp>("case_t0",case_t0);
  c.create<TimeStamp>("stop_timestamp",stop);

  if(params.isSublist("Profiling Output")) {
    c.create<ProfilingArchive>("archive",comm,case_t0,run_t0,params.sublist("Profiling Output"));
  } else {
    ekat::ParameterList profiling_output_list("Profiling Output");
    c.create<ProfilingArchive>("archive",comm,case_t0,run_t0,profiling_output_list);
  }

  if (comm.am_i_root()) {
    printf(" [CLDERA] Initializing profiling context '%s' ... done!\n",context_name);
  }
  c.timing().stop_timer(timer_name);
}

void cldera_clean_up_c ()
{
  auto& c = get_curr_context();

  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  // Need a copy for the last print statement, when the context won't be available anymore
  const auto cname = c.name();

  // Store this, since cleaning up the ProfileContext will reset
  // the comm to MPI_COMM_SELF.
  auto am_i_root = c.get_comm().am_i_root();
  if (am_i_root) {
    printf(" [CLDERA] Shutting down profiling context '%s' ...\n", cname.c_str());
  }

  auto& params = c.get_params();
  if(params.isSublist("Pathway")) {
    const auto& history_filename = params.get<std::string>("pathway_history_file","cldera_pathway_history.yaml");
    auto& pathway = c.get<std::shared_ptr<cldera::Pathway>>("pathway");
    pathway->dump_test_history_to_yaml(history_filename);
  }

  // Clean up and remove context from the session
  // If this was the last stored context, this call will also clean up the session
  get_session().clean_up_current_context();

  if (am_i_root) {
    printf(" [CLDERA] Shutting down profiling context '%s'... done!\n",cname.c_str());
  }
}

void cldera_add_field_c (const char*& name,
                         const int    rank,
                         const int*   dims,
                         const char** dimnames,
                         const bool   is_view,
                         const char*& data_type)
{
  cldera_add_partitioned_field_c(name,rank,dims,dimnames,1,0,-1,is_view,data_type);
}

void cldera_add_partitioned_field_c (
    const char*&  name,
    const int     rank,
    const int*    dims,
    const char**  dimnames,
    const int     num_parts,
    const int     part_dim,
    const int     part_dim_alloc_size,
    const bool    is_view,
    const char*&  dtype)
{
  auto& c = get_curr_context();

  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  auto& ts = c.timing();
  ts.start_timer(c.name() + "::add_field");

  EKAT_REQUIRE_MSG (rank>=0 && rank<=4,
      "Error! Unsupported field rank (" + std::to_string(rank) + "\n");
  EKAT_REQUIRE_MSG (num_parts>=1,
      "Error! Invalid number of partitions (" + std::to_string(num_parts) + "\n");

  // Copy input raw pointer to vector
  std::vector<int> d(rank);
  std::vector<std::string> dn(rank);
  for (int i=0; i<rank; ++i) {
    EKAT_REQUIRE_MSG (dims[i]>=0,
        "Error! Invalid field extent.\n"
        "   - Field name:" +  std::string(name) + "\n";
        "   - Dimension: " +  std::to_string(i) + "\n"
        "   - Extent:    " +  std::to_string(dims[i]) + "\n");

    d[i] = dims[i];
    dn[i] = dimnames[i];
  }
  FieldLayout fl(d,dn);

  // Set data in the archive structure
  auto& archive = c.get<ProfilingArchive>("archive");
  const auto access = is_view ? DataAccess::View : DataAccess::Copy;
  Field f(name,fl,num_parts,part_dim,access,str2data_type(dtype),part_dim_alloc_size);
  archive.add_field(f);
  ts.stop_timer(c.name() + "::add_field");
}

void cldera_set_field_part_extent_c (
    const char*& name,
    const int   part,
    const int   part_extent)
{
  auto& c = get_curr_context();

  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  auto& ts = c.timing();
  ts.start_timer(c.name() + "::set_field_size");

  auto& archive = c.get<ProfilingArchive>("archive");

  archive.get_field(name).set_part_extent (part,part_extent);
  ts.stop_timer(c.name() + "::set_field_size");
}

void cldera_set_field_part_data_c (
    const char*& name,
    const int   part,
    const void*& data_in,
    const char*& dtype_in)
{
  auto& c = get_curr_context();

  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  auto& ts = c.timing();
  ts.start_timer(c.name() + "::set_field_data");

  auto& archive = c.get<ProfilingArchive>("archive");

  auto& f = archive.get_field(name);
  std::string dtype = dtype_in;
  if (dtype=="real") {
    const Real* data = reinterpret_cast<const Real*>(data_in);
    if (f.data_access()==DataAccess::View) {
      f.set_part_data<Real> (part,data);
    } else {
      f.copy_part_data<Real> (part,data);
    }
  } else if (dtype=="int") {
    const int* data = reinterpret_cast<const int*>(data_in);
    if (f.data_access()==DataAccess::View) {
      f.set_part_data<int> (part,data);
    } else {
      f.copy_part_data<int> (part,data);
    }
  } else {
    EKAT_ERROR_MSG ("Invalid/unsupported data type: " + dtype + "\n");
  }
  ts.stop_timer(c.name() + "::set_field_data");
}

void cldera_commit_field_c (const char*& name)
{
  auto& c = get_curr_context();

  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  auto& ts = c.timing();
  ts.start_timer(c.name() + "::commit_fields");

  auto& archive = c.get<ProfilingArchive>("archive");

  archive.get_field(name).commit();
  ts.stop_timer(c.name() + "::commit_fields");
}

void cldera_commit_all_fields_c ()
{
  auto& c = get_curr_context();

  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  auto& ts = c.timing();
  auto& archive = c.get<ProfilingArchive>("archive");

  ts.start_timer(c.name() + "::commit_fields");
  archive.commit_all_fields();
  ts.stop_timer(c.name() + "::commit_fields");

  ts.start_timer(c.name() + "::create_stats");
  auto& params = c.get_params();
  using stat_ptr_t = std::shared_ptr<FieldStat>;
  using requests_t = std::map<std::string,std::vector<stat_ptr_t>>;
  using vos_t = std::vector<std::string>;

  auto& factory = StatFactory::instance();
  register_stats();
  auto& requests = c.create<requests_t>("requests");
  const auto& fnames = params.get<vos_t>("Fields To Track");
  for (const auto& fname : fnames) {
    auto& req_pl = params.sublist(fname);
    auto& req_stats = requests[fname];
    const auto& f = archive.get_field(fname);
    for (auto stat_name : req_pl.get<vos_t>("Compute Stats")) {
      auto& stat_pl = req_pl.sublist(stat_name);
      const auto& stat_type = stat_pl.get<std::string>("type",stat_name);
      auto stat = factory.create(stat_type,c.get_comm(),stat_pl);
      stat->set_field(f);
      std::map<std::string,Field> aux_fields;
      for (const auto& fn : stat->get_aux_fields_names()) {
        // It may seem that we should let archive error out here,
        // if fn is not found. However, there'c a good case for not doing it.
        // Namely, all masked integral stats need a mask field. When the first
        // instance is created, we have no mask in the archive, so the stat
        // will load from nc file. After that, we'll put the stat in the archive
        // so that we can pass it to the other instances using the same mask
        if (archive.has_field(fn)) {
          aux_fields[fn] = archive.get_field(fn);
        }
      }
      stat->set_aux_fields(aux_fields);
      stat->create_stat_field();
      req_stats.push_back(stat);

      // Add all fields computed by the stat to the archive
      if (not archive.has_field(stat->get_stat_field().name())) {
        archive.add_field(stat->get_stat_field());
      }
      for (const auto& it : stat->get_aux_fields()) {
        if (not archive.has_field(it.first)) {
          archive.add_field(it.second);
        }
      }
    }
  }
  ts.stop_timer(c.name() + "::create_stats");
}

void cldera_compute_stats_c (const int ymd, const int tod)
{
  auto& c = get_curr_context();
  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  cldera::TimeStamp time = {ymd, tod};
  if (time==c.get<TimeStamp>("run_t0")) {
    // E3SM runs a bit of its timestep during init, to bootstrap
    // some surface fluxes. We are not interested in the stats at
    // that time.
    return;
  }
  static std::map<std::string,int> num_calls_map;
  auto& num_calls = num_calls_map[c.name()];

  const auto& comm = c.get_comm();
  auto& params = c.get_params();

  // There is some rank that enters this call after the others,
  // making the relative timing of internal cldera funcs harder.
  // We can use the following to add a barrier upon entrance in
  // this routine, so that all ranks will start together
  if (params.get<bool>("Add Compute Stats Barrier",false)) {
    comm.barrier();
  }

  if (comm.am_i_root()) {
    printf(" [CLDERA] Computing stats for context '%s'...\n",c.name().c_str());
    printf(" [CLDERA]   time: %s...\n",time.to_string().c_str());
  }

  auto& ts = c.timing();
  ts.start_timer(c.name() + "::compute_stats");

  using stat_ptr_t = std::shared_ptr<FieldStat>;
  using requests_t = std::map<std::string,std::vector<stat_ptr_t>>;
  auto& requests = c.get<requests_t>("requests");

  auto& archive = c.get<ProfilingArchive>("archive");

  std::map<std::string, std::shared_ptr<const cldera::Field> > fields;
  for (const auto& it : requests) {
    const auto& fname = it.first;
    const auto& stats = it.second;
    const auto& f = archive.get_field(fname);

    // Store field map, so we can create the pathway object later
    if(!c.has_data("pathway")) {
      ts.start_timer("profiling:create_stat_field");
      fields[fname] = std::make_shared<const cldera::Field>(f);
      ts.stop_timer("profiling:create_stat_field");
    }
    for (auto stat : stats) {
      archive.append_stat(fname,stat->name(),stat->compute(time));
    }
  }

  archive.update_time(time);
  ts.stop_timer(c.name() + "::compute_stats");

  if (comm.am_i_root()) {
    printf(" [CLDERA] Computing stats for context '%s'...done!\n",c.name().c_str());
  }

  // If pathway is enabled, do it
  if (params.isSublist("Pathway")) {
    c.timing().start_timer(c.name() + "::run_pathway_tests");
    // this solves the issue of fields not being initialized
    if(not c.has_data("pathway")) {
      auto& archive = c.get<ProfilingArchive>("archive");
      // It's the first time this is called, so create the pathway
      cldera::PathwayFactory pathway_factory(params, comm, false); // TODO: make pathway verbosity a yaml option
      c.create<std::shared_ptr<Pathway>>("pathway",pathway_factory.build_pathway(archive));
    }
    auto pathway = c.get<std::shared_ptr<cldera::Pathway>>("pathway");
    pathway->run_pathway_tests(comm, time);
    c.timing().stop_timer(c.name() + "::run_pathway_tests");
  }

  const int timings_flush_freq = params.get("Timings Flush Freq",0);
  if (ts.is_active() and timings_flush_freq>0 and num_calls%timings_flush_freq==0) {
    const auto timings_fname = params.get<std::string>("Timing Filename","");
    const auto& timing = c.timing();
    std::ofstream timings_file;
    std::stringstream blackhole;
    if (comm.am_i_root()) {
      timings_file.open(timings_fname);
    }
    std::ostream& ofile = timings_file;
    std::ostream& onull = blackhole;

    std::ostream& out = comm.am_i_root() ? ofile : onull;
    timing.dump(out,comm);
  }
  ++num_calls;
}

} // namespace cldera

} // extern "C"
