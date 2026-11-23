#include "cldera_profiling_interface.hpp"

#include "cldera_profiling_context.hpp"
#include "cldera_profiling_session.hpp"
#include "cldera_profiling_archive.hpp"
#include "cldera_pathway_factory.hpp"
#include "stats/cldera_register_stats.hpp"
#include "stats/cldera_field_single_rank.hpp"

#include "timing/cldera_timing_session.hpp"

#include "cldera_sai_writer.hpp"

#ifdef CLDERA_CONTROLLER
  #include "cldera_controller.hpp"
  #include "cldera_controller_python.hpp"
  #include "cldera_annual_mean.hpp"
#endif

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

#ifdef CLDERA_CONTROLLER
  if (params.isParameter("Controller")) {
    // BMW: pass the case start, so that the controller yaml may use 'Start Years After Case Start'
    auto& ctrl = c.create<Controller>("controller", params, case_t0);

    if (comm.am_i_root()) {
      printf(" [CLDERA] Controller '%s': active from %s (%s) to %s\n",
             ctrl.get_filename().c_str(),
             TimeStamp(ctrl.get_start_ymd(),ctrl.get_start_tod()).to_string().c_str(),
             ctrl.get_start_mode().c_str(),
             TimeStamp(ctrl.get_end_ymd(),ctrl.get_end_tod()).to_string().c_str());
    }

    // BMW (G6-1.5K-HiLLA): annual means of the feedback stats, for annual feedback
    //      controllers. Accumulated at every compute_stats call; state is carried across
    //      restarts with a small text file (see cldera_annual_mean.hpp).
    const auto fb_stats = ctrl.get_feedback_stats();
    if (not fb_stats.empty()) {
      auto& p = c.get_params();
      const auto prefix = p.get<std::string>("Annual Mean Restart File Prefix","cldera_annual_mean.rst");
      const bool required = p.get<bool>("Annual Mean Restart Required",true);
      c.create<AnnualMeanAccumulator>("annual_mean",comm,fb_stats,prefix,required,case_t0,run_t0);
    }

    // BMW: if the context yaml has an 'SAI Writer' sublist, create the C++ writer that
    //      applies the controller output to the (writable) EAM 'sai' forcing field.
    //      Sites (lat/lon/alt) come from the controller yaml 'Injections', so they are
    //      defined in ONE place only.
    if (params.isSublist("SAI Writer")) {
      std::vector<SAIWriter::Site> sites;
      const auto names = ctrl.get_inj_names();
      const auto lats  = ctrl.get_inj_lats();
      const auto lons  = ctrl.get_inj_lons();
      const auto alts  = ctrl.get_inj_alts();
      for (size_t i=0; i<names.size(); ++i) {
        // Controller yaml 'alt' is in km
        sites.push_back(SAIWriter::Site{names[i],lats[i],lons[i],alts[i]*1000.0});
      }
      c.create<SAIWriter>("sai_writer",comm,params.sublist("SAI Writer"),sites);
    }
  }
#endif

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

  // Handle Python cleanup here to avoid shutdown ordering errors
#ifdef CLDERA_CONTROLLER
  // BMW: save the annual-mean accumulator state for the next run (named by this run's stop
  //      time, which is the next run's start time), and close the SAI writer's mass budget.
  if (c.has_data("annual_mean")) {
    c.get<AnnualMeanAccumulator>("annual_mean").write_restart(c.get<TimeStamp>("stop_timestamp"));
  }
  if (c.has_data("sai_writer")) {
    c.get<SAIWriter>("sai_writer").finalize();
  }
  controller_finalize_python();
#endif

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

void cldera_set_field_part_data_nonconst_c (
    const char*& name,
    const int   part,
    void*& data_in,
    const char*& dtype_in)
{
  // Like cldera_set_field_part_data_c, but hands cldera a NON-CONST pointer to
  // host-app memory, so that cldera may write into it. Only valid for View fields
  // registered as writable (see cldera_add_partitioned_field_c, read_only=false).
  auto& c = get_curr_context();

  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  auto& ts = c.timing();
  ts.start_timer(c.name() + "::set_field_data");

  auto& archive = c.get<ProfilingArchive>("archive");

  auto& f = archive.get_field(name);
  EKAT_REQUIRE_MSG (f.data_access()==DataAccess::View and not f.is_read_only(),
      "[cldera_set_field_part_data_nonconst_c] Error! Field must be a writable View field.\n"
      " - field name: " + std::string(name) + "\n");
  std::string dtype = dtype_in;
  if (dtype=="real") {
    f.set_part_data<Real> (part,reinterpret_cast<Real*>(data_in));
  } else if (dtype=="int") {
    f.set_part_data<int> (part,reinterpret_cast<int*>(data_in));
  } else {
    EKAT_ERROR_MSG ("Invalid/unsupported data type: " + dtype + "\n");
  }
  ts.stop_timer(c.name() + "::set_field_data");
}

void cldera_set_field_writable_c (const char*& name)
{
  auto& c = get_curr_context();
  if (not c.inited()) { return; }

  auto& archive = c.get<ProfilingArchive>("archive");
  archive.get_field(name).set_read_only(false);
}

void cldera_set_sai_storage_bottom_up_c (const bool bottom_up)
{
  auto& c = get_curr_context();
  if (not c.inited()) { return; }

  // No-op if no SAI writer was requested in the yaml
  if (c.has_data("sai_writer")) {
    c.get<SAIWriter>("sai_writer").set_storage_bottom_up(bottom_up);
    if (c.get_comm().am_i_root()) {
      printf(" [CLDERA] SAIWriter: target field storage is %s\n",
             bottom_up ? "BOTTOM-UP (EAM alt_data, will be flipped by extfrc_set)" : "TOP-DOWN (model order)");
    }
  }
}

int cldera_has_sai_writer_c ()
{
  auto& c = get_curr_context();
  if (not c.inited()) { return 0; }
  return c.has_data("sai_writer") ? 1 : 0;
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
        // if fn is not found. However, there's a good case for not doing it.
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

  // Cache static SingleRank stats for gids/area once, so all ranks participate
#ifdef CLDERA_CONTROLLER
  if (!c.has_data("controller_static_singlerank")) {
    auto& cached = c.create<std::vector<std::shared_ptr<FieldSingleRank>>>("controller_static_singlerank");
    if (archive.has_field("col_gids")) {
      ekat::ParameterList pl("col_gids");
      auto sr = std::make_shared<FieldSingleRank>(comm, pl);
      sr->set_field(archive.get_field("col_gids"));
      sr->create_stat_field();
      sr->compute(time);
      cached.push_back(sr);
    }
    if (archive.has_field("area")) {
      ekat::ParameterList pl("area");
      auto sr = std::make_shared<FieldSingleRank>(comm, pl);
      sr->set_field(archive.get_field("area"));
      sr->create_stat_field();
      sr->compute(time);
      cached.push_back(sr);
    }
  }
#endif

  for (const auto& it : requests) {
    const auto& fname = it.first;
    const auto& stats = it.second;
    const auto& f = archive.get_field(fname);

    for (auto& stat : stats) {
      auto stat_field = stat->compute(time);
      // Do not write singlerank stats to disk (avoid heavy I/O and rank-dim issues).
      if (stat->type() == "singlerank") {
        continue;
      }
      archive.update_stat(fname,stat->name(),stat_field);
    }
  }

  archive.end_timestep(time);
  ts.stop_timer(c.name() + "::compute_stats");

  if (comm.am_i_root()) {
    printf(" [CLDERA] Computing stats for context '%s'...done!\n",c.name().c_str());
  }

  // BMW: controller + SAI writer. The controller runs (at most) once per control
  //      interval, using the stats just computed. The SAI writer runs EVERY call,
  //      because EAM re-fills the 'sai' forcing array from file (zeros) at every
  //      step in interpolate_trcdata, before tphysbc; this call happens after
  //      tphysbc and before tphysac/chemistry, which is where extfrc_set reads it.
#ifdef CLDERA_CONTROLLER
  // BMW (G6-1.5K-HiLLA): add this step's feedback stats to the annual means BEFORE the
  //      controller runs. The Jan 1 00:00 sample opens the new year, which closes the
  //      previous one, so the controller triggered at Jan 1 00:00 sees the full
  //      previous year. Global stats are all-reduced, so root has the right values.
  if (c.has_data("annual_mean")) {
    auto& am = c.get<AnnualMeanAccumulator>("annual_mean");
    std::vector<double> vals;
    for (const auto& sname : am.stat_names()) {
      bool found = false;
      for (auto& it : requests) {
        for (auto& stat : it.second) {
          if (stat->name()!=sname) continue;
          const auto& sf = stat->get_stat_field();
          EKAT_REQUIRE_MSG (sf.committed() and sf.data_type()==RealType and sf.nparts()==1
                            and sf.layout().size()==1,
              "[CLDERA] Error! Controller feedback stat '" + sname + "' must be a committed\n"
              "  scalar Real stat (e.g. a masked_integral over the global mask).\n");
          vals.push_back(sf.view<const Real>()(0));
          found = true;
        }
      }
      EKAT_REQUIRE_MSG (found,
          "[CLDERA] Error! Controller feedback stat '" + sname + "' not found in the requests.\n");
    }
    am.add(time,vals);
  }
  cldera_compute_controller_c(ymd,tod);
#endif
  if (c.has_data("sai_writer")) {
    ts.start_timer(c.name() + "::sai_writer");
    c.get<SAIWriter>("sai_writer").apply(archive,time);
    ts.stop_timer(c.name() + "::sai_writer");
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


#ifdef CLDERA_CONTROLLER
void cldera_compute_controller_c (const int ymd, const int tod)
{
  auto& c = get_curr_context();
  // If input file was not provided, cldera does nothing
  if (not c.inited()) { return; }

  const auto& comm = c.get_comm();
  auto& params = c.get_params();

  // Case is not controlled
  if (!params.isParameter("Controller")) { return; }

  auto& ctrl = c.get<Controller>("controller");

  // BMW: trigger at the start of each control interval (tod==0 on the interval's
  //      start day), AND on the very first call of each run. The latter is needed
  //      because (a) the first stats call of a startup run is at tod>0 (E3SM
  //      #5904 workaround), and (b) on restarts the controller history is not
  //      saved, so we must recompute the current rates.
  //      Interval starts trigger the Python controller also OUTSIDE the controller
  //      active window [Start, End]: there it returns zeros, so that injection is off
  //      before the start and stops at the end (otherwise the last rates would be kept).
  const bool first_call = not c.has_data("ctrl_mass_inj");
  if (not first_call and !ctrl.is_interval_start(ymd, tod)) { return; }
  // Guard against double triggering on the same day (Controller::operator() errors out on that).
  // Can happen if the first call of the run is at tod==0 of an interval start day.
  if (ctrl.has_ctrl_hist_at(ymd)) { return; }

  cldera::TimeStamp time = {ymd, tod};
  if (comm.am_i_root()) {
    printf(" [CLDERA] Computing controller for context '%s'%s...\n",c.name().c_str(),
           first_call ? " (first call of this run)" : "");
    printf(" [CLDERA]   time: %s...\n",time.to_string().c_str());
  }

  auto& ts = c.timing();
  ts.start_timer(c.name() + "::compute_controller");

  // Call controller (Python runs on root only; other ranks get an empty vector)
  ctrl.set_context(&c);
  ctrl.set_call_info(first_call, c.get<TimeStamp>("run_t0"));
  std::vector<double> mass_inj = ctrl(ymd, tod);

  // BMW: broadcast the result to all ranks. Every rank knows the number of sites
  //      from the yaml, so the size is known everywhere. Check the size on root first,
  //      and broadcast the check, so that all ranks abort consistently.
  const int ninj = ctrl.get_inj_names().size();
  int ok = comm.am_i_root() ? (static_cast<int>(mass_inj.size())==ninj ? 1 : 0) : 1;
  comm.broadcast(&ok,1,comm.root_rank());
  EKAT_REQUIRE_MSG (ok==1,
      "[cldera_compute_controller_c] Error! Python controller returned the wrong number of values.\n"
      " - expected (number of 'Injections'): " + std::to_string(ninj) + "\n");
  mass_inj.resize(ninj);
  comm.broadcast(mass_inj.data(),ninj,comm.root_rank());
  // Record the broadcast value in the history on ALL ranks (operator() only has it on root),
  // so that has_ctrl_hist_at behaves identically everywhere.
  ctrl.set_ctrl_hist_at(ymd,mass_inj);

  // Store result, so the SAI writer (and anyone else) can use it until the next trigger
  if (first_call) {
    c.create<std::vector<double>>("ctrl_mass_inj",mass_inj);
  } else {
    c.get<std::vector<double>>("ctrl_mass_inj") = mass_inj;
  }
  if (c.has_data("sai_writer")) {
    // Interval name, for the per-interval mass budget in the SAI writer
    std::string iname;
    const int mmdd = ymd % 10000;
    auto doy = [](const int md) { return day_of_year_365(md/100, md%100); };
    for (const auto& idef : ctrl.get_intervals()) {
      const int s = doy(idef.start_mmdd), e = doy(idef.end_mmdd), d = doy(mmdd);
      const bool in = e>=s ? (d>=s and d<=e) : (d>=s or d<=e);
      if (in) iname = idef.name;
    }
    c.get<SAIWriter>("sai_writer").set_rates(mass_inj,time,iname);
  }

  ts.stop_timer(c.name() + "::compute_controller");
}

void cldera_get_controller_output_c (double* vals, int* nvals, int* has_value)
{
  auto& c = get_curr_context();
  *has_value = 0;
  if (not c.inited() or not c.has_data("ctrl_mass_inj")) { *nvals = 0; return; }
  const auto& v = c.get<std::vector<double>>("ctrl_mass_inj");
  EKAT_REQUIRE_MSG (*nvals >= static_cast<int>(v.size()),
      "[cldera_get_controller_output_c] Error! Output buffer too small.\n");
  std::copy(v.begin(),v.end(),vals);
  *nvals = v.size();
  *has_value = 1;
}
#else
// Stubs, so that the host app can link even without the controller enabled
void cldera_compute_controller_c (const int /* ymd */, const int /* tod */) {}
void cldera_get_controller_output_c (double* /* vals */, int* nvals, int* has_value)
{
  *nvals = 0;
  *has_value = 0;
}
#endif

} // namespace cldera

} // extern "C"
