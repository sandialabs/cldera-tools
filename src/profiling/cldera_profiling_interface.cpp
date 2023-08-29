#include "cldera_profiling_interface.hpp"

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

void cldera_init_c (const MPI_Fint fcomm,
                    const int case_t0_ymd, const int case_t0_tod,
                    const int run_t0_ymd, const int run_t0_tod,
                    const int stop_ymd, const int stop_tod)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::init");

  // Convert F90 comm to C comm, create ekat::Comm, and init session
  MPI_Comm mpiComm = MPI_Comm_f2c(fcomm);
  EKAT_REQUIRE_MSG (mpiComm!=nullptr, "Error! Input fortran comm seems invalid.\n");

  ekat::Comm comm(mpiComm);

  ekat::initialize_ekat_session (/* argc = */ 0,
                                 /* argv = */ nullptr,
                                 /* print config = */ comm.am_i_root());

  //TODO: make the filename configurable
  std::string filename = "./cldera_profiling_config.yaml";

  if (not std::ifstream(filename).good()) {
    if (comm.am_i_root()) {
      printf(" [CLDERA] WARNING: could not open './cldera_profiling_config.yaml'.\n"
             "   -> Profiling will do nothing.\n");
    }
    ts.stop_timer("profiling::init");
    return;
  }

  if (comm.am_i_root()) {
    printf(" [CLDERA] Initializing profiling session ...\n");
  }

  auto& s = ProfilingSession::instance();

  s.init(comm);

  using stat_ptr_t = std::shared_ptr<FieldStat>;
  using requests_t = std::map<std::string,std::vector<stat_ptr_t>>;
  using vos_t = std::vector<std::string>;

  auto& params = s.get_params() = ekat::parse_yaml_file(filename);
  TimeStamp case_t0 (case_t0_ymd,case_t0_tod);
  TimeStamp run_t0 (run_t0_ymd,run_t0_tod);
  TimeStamp stop (stop_ymd,stop_tod);
  s.create<TimeStamp>("run_t0",run_t0);
  s.create<TimeStamp>("case_t0",case_t0);
  s.create<TimeStamp>("stop_timestamp",stop);

  if(params.isSublist("Profiling Output")) {
    s.create<ProfilingArchive>("archive",comm,case_t0,run_t0,params.sublist("Profiling Output"));
  } else {
    ekat::ParameterList profiling_output_list("Profiling Output");
    s.create<ProfilingArchive>("archive",comm,case_t0,run_t0,profiling_output_list);
  }

  s.create<bool>("doPathway",params.isSublist("Pathway"));

  if (comm.am_i_root()) {
    printf(" [CLDERA] Initializing profiling session ... done!\n");
  }
  ts.stop_timer("profiling::init");
}

void cldera_clean_up_c ()
{
  auto& s = ProfilingSession::instance();

  // If input file was not provided, cldera does nothing
  if (not s.inited()) { return; }

  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::clean_up");

  // Store this, since cleaning up the ProfileSession will reset
  // the comm to MPI_COMM_SELF.
  auto am_i_root = s.get_comm().am_i_root();
  if (am_i_root) {
    printf(" [CLDERA] Shutting down profiling session ...\n");
  }

  if(s.get<bool>("doPathway")) {
    std::string history_filename = "cldera_pathway_history.yaml";
    auto& pathway = s.get<std::shared_ptr<cldera::Pathway>>("pathway");
    pathway->dump_test_history_to_yaml(history_filename);
  }

  // Get a copy of timing filename and comm *before* cleaning up the
  // session, since we need them to dump the timing stats
  auto& params = s.get_params();
  const auto timings_fname = params.get<std::string>("Timing Filename","");
  const auto comm = s.get_comm();
  s.clean_up();
  ts.stop_timer("profiling::clean_up");

  if (timings_fname!="") {
    const auto& timings = timing::TimingSession::instance();
    std::ofstream timings_file;
    std::stringstream blackhole;
    if (am_i_root) {
      timings_file.open(timings_fname);
    }
    std::ostream& ofile = timings_file;
    std::ostream& onull = blackhole;

    std::ostream& out = am_i_root ? ofile : onull;
    timings.dump(out,comm);
  }

  if (am_i_root) {
    printf(" [CLDERA] Shutting down profiling session ... done!\n");
  }
  ekat::finalize_ekat_session ();
}

void cldera_add_field_c (const char*& name,
                         const int    rank,
                         const int*   dims,
                         const char** dimnames,
                         const bool   is_view,
                         const char*& data_type)
{
  cldera_add_partitioned_field_c(name,rank,dims,dimnames,1,0,is_view,data_type);
}

void cldera_add_partitioned_field_c (
    const char*&  name,
    const int     rank,
    const int*    dims,
    const char**  dimnames,
    const int     num_parts,
    const int     part_dim,
    const bool    is_view,
    const char*&  dtype)
{
  auto& s = ProfilingSession::instance();

  // If input file was not provided, cldera does nothing
  if (not s.inited()) { return; }

  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::add_field");

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
  auto& archive = s.get<ProfilingArchive>("archive");
  const auto access = is_view ? DataAccess::View : DataAccess::Copy;
  archive.add_field(Field(name,fl,num_parts,part_dim,access,str2data_type(dtype)));
  ts.stop_timer("profiling::add_field");
}

void cldera_set_field_part_size_c (
    const char*& name,
    const int   part,
    const int   part_size)
{
  auto& s = ProfilingSession::instance();

  // If input file was not provided, cldera does nothing
  if (not s.inited()) { return; }

  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::set_field_size");

  auto& archive = s.get<ProfilingArchive>("archive");

  archive.get_field(name).set_part_size (part,part_size);
  ts.stop_timer("profiling::set_field_size");
}

void cldera_set_field_part_data_c (
    const char*& name,
    const int   part,
    const void*& data_in,
    const char*& dtype_in)
{
  auto& s = ProfilingSession::instance();

  // If input file was not provided, cldera does nothing
  if (not s.inited()) { return; }

  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::set_field_data");

  auto& archive = s.get<ProfilingArchive>("archive");

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
  ts.stop_timer("profiling::set_field_data");
}

void cldera_commit_field_c (const char*& name)
{
  auto& s = ProfilingSession::instance();

  // If input file was not provided, cldera does nothing
  if (not s.inited()) { return; }

  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::commit_fields");

  auto& archive = s.get<ProfilingArchive>("archive");

  archive.get_field(name).commit();
  ts.stop_timer("profiling::commit_fields");
}

void cldera_commit_all_fields_c ()
{
  auto& s = ProfilingSession::instance();

  // If input file was not provided, cldera does nothing
  if (not s.inited()) { return; }

  auto& ts = timing::TimingSession::instance();
  auto& archive = s.get<ProfilingArchive>("archive");

  ts.start_timer("profiling::commit_fields");
  archive.commit_all_fields();
  ts.stop_timer("profiling::commit_fields");

  ts.start_timer("profiling::create_stats");
  auto& params = s.get_params();
  using stat_ptr_t = std::shared_ptr<FieldStat>;
  using requests_t = std::map<std::string,std::vector<stat_ptr_t>>;
  using vos_t = std::vector<std::string>;

  auto& factory = StatFactory::instance();
  register_stats();
  auto& requests = s.create<requests_t>("requests");
  const auto& fnames = params.get<vos_t>("Fields To Track");
  for (const auto& fname : fnames) {
    auto& req_pl = params.sublist(fname);
    auto& req_stats = requests[fname];
    const auto& f = archive.get_field(fname);
    for (auto stat_name : req_pl.get<vos_t>("Compute Stats")) {
      auto& stat_pl = req_pl.sublist(stat_name);
      const auto& stat_type = stat_pl.get<std::string>("type",stat_name);
      auto stat = factory.create(stat_type,s.get_comm(),stat_pl);
      stat->set_field(f);
      std::map<std::string,Field> aux_fields;
      for (const auto& fn : stat->get_aux_fields_names()) {
        aux_fields[fn] = archive.get_field(fn);
      }
      stat->set_aux_fields(aux_fields);
      stat->create_stat_field();
      req_stats.push_back(stat);
    }
  }
  ts.stop_timer("profiling::create_stats");
}

void cldera_compute_stats_c (const int ymd, const int tod)
{
  auto& s = ProfilingSession::instance();
  // If input file was not provided, cldera does nothing
  if (not s.inited()) { return; }

  cldera::TimeStamp time = {ymd, tod};
  if (time==s.get<TimeStamp>("run_t0")) {
    // E3SM runs a bit of its timestep during init, to bootstrap
    // some surface fluxes. We are not interested in the stats at
    // that time.
    return;
  }
  static int num_calls = 0;

  const auto& comm = s.get_comm();
  auto& params = s.get_params();

  // There is some rank that enters this call after the others,
  // making the relative timing of internal cldera funcs harder.
  // We can use the following to add a barrier upon entrance in
  // this routine, so that all ranks will start together
  if (params.get<bool>("Add Compute Stats Barrier",false)) {
    comm.barrier();
  }

  if (comm.am_i_root()) {
    printf(" [CLDERA] Computing stats at time %s...\n",time.to_string().c_str());
  }

  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::compute_stats");

  using stat_ptr_t = std::shared_ptr<FieldStat>;
  using requests_t = std::map<std::string,std::vector<stat_ptr_t>>;
  auto& requests = s.get<requests_t>("requests");

  auto& archive = s.get<ProfilingArchive>("archive");

  std::map<std::string, std::shared_ptr<const cldera::Field> > fields;
  for (const auto& it : requests) {
    const auto& fname = it.first;
    const auto& stats = it.second;
    const auto& f = archive.get_field(fname);

    // Store field map, so we can create the pathway object later
    if(!s.has_data("pathway")) {
      ts.start_timer("profiling:create_stat_field");
      fields[fname] = std::make_shared<const cldera::Field>(f);
      ts.stop_timer("profiling:create_stat_field");
    }
    for (auto stat : stats) {
      archive.append_stat(fname,stat->name(),stat->compute(time));
    }
  }

  archive.update_time(time);
  ts.stop_timer("profiling::compute_stats");

  if (comm.am_i_root()) {
    printf(" [CLDERA] Computing stats ... done!\n");
  }

  if(s.get<bool>("doPathway")) {
    ts.start_timer("profiling::run_pathway_tests");
    // this solves the issue of fields not being initialized
    if(!s.has_data("pathway")) {
      // create the pathway then throw it in the ProfilingSession
      std::string filename = "./cldera_profiling_config.yaml";
      cldera::PathwayFactory pathway_factory(filename, fields, comm, false); // TODO: make pathway verbosity a yaml option
      auto& pathway = s.create_or_get<std::shared_ptr<cldera::Pathway>>("pathway",pathway_factory.build_pathway());
      pathway->run_pathway_tests(comm, time);
    } else {
      // otherwise, grab it
      auto& pathway = s.get<std::shared_ptr<cldera::Pathway>>("pathway");
      pathway->run_pathway_tests(comm, time);
    }
    ts.stop_timer("profiling::run_pathway_tests");
  }

  const auto timings_fname = params.get<std::string>("Timing Filename","");
  const int timings_flush_freq = params.get("Timings Flush Freq",0);
  if (timings_fname!="" && timings_flush_freq>0 && num_calls%timings_flush_freq==0) {
    const auto& timings = timing::TimingSession::instance();
    std::ofstream timings_file;
    std::stringstream blackhole;
    if (comm.am_i_root()) {
      timings_file.open(timings_fname);
    }
    std::ostream& ofile = timings_file;
    std::ostream& onull = blackhole;

    std::ostream& out = comm.am_i_root() ? ofile : onull;
    timings.dump(out,comm);
  }
  ++num_calls;
}

} // namespace cldera

} // extern "C"
