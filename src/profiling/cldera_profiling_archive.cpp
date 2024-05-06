#include "cldera_profiling_archive.hpp"
#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"
#include "timing/cldera_timing_session.hpp"

#include <ekat/io/ekat_yaml.hpp>
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <numeric>
#include <fstream>

namespace cldera {

ProfilingArchive::
ProfilingArchive(const ekat::Comm& comm,
                 const TimeStamp& case_t0,
                 const TimeStamp& run_t0,
                 const ekat::ParameterList& params)
 : m_comm (comm)
 , m_params (params)
 , m_case_t0 (case_t0)
{
  using intvec_t = std::vector<int>;
  if (m_params.get<bool>("Enable Output",true)) {
    const auto& prefix = m_params.get<std::string>("filename_prefix","cldera_stats");
    std::string filename;
    io::pnetcdf::IOMode mode = io::pnetcdf::IOMode::Invalid;

    const auto& time_avg_sizes = m_params.get<intvec_t>("time_averaging_window_sizes",intvec_t(1,1));
    m_num_streams = time_avg_sizes.size();
    for (auto s : time_avg_sizes) {

      auto suffix = s==1 ? std::string(".INSTANT") : ".AVERAGE.nsteps_x" + std::to_string(s);
      if (m_case_t0==run_t0) {
        filename = prefix + "-" + m_case_t0.to_string();
        mode = io::pnetcdf::IOMode::Write;
      } else {
        filename = prefix + "-" + m_case_t0.to_string();
        mode = io::pnetcdf::IOMode::Append;
        if (not std::ifstream(filename).good() or m_params.get("force_new_file",true)) {
          filename = prefix + "-" + run_t0.to_string();
          mode = io::pnetcdf::IOMode::Write;
        }
      }
      m_output_files.push_back(open_file (filename+".nc", m_comm, mode));

      // Init the beg/end vector of each averaging window
      m_time_avg_beg.resize(m_num_streams,run_t0);
      m_time_avg_end.resize(m_num_streams);

      m_time_avg_window_size.push_back(s);
      m_time_avg_curr_count.push_back(0);
    }
    m_fields_stats.resize(m_num_streams);
  }
}

ProfilingArchive::
~ProfilingArchive()
{
  for (int i=0; i<m_num_streams; ++i) {
    auto& f = m_output_files[i];
    if (f) {
      io::pnetcdf::close_file(*f);
    }
  }
}

void ProfilingArchive::
setup_output_file (const int istream)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::setup_output_file");

  if (m_comm.am_i_root()) {
    printf(" [CLDERA] setting up output file ...");
  }

  auto& file = *m_output_files[istream];

  // Add time dimension
  io::pnetcdf::add_time (file,"double");

  if (m_time_avg_window_size[istream]>1) {
    io::pnetcdf::add_dim(file, "dim2", 2);
    // We save the start/end of each averaging window
    io::pnetcdf::add_var (file,
                          "time_bounds",
                          io::pnetcdf::get_io_dtype_name<Real>(),
                          {"dim2"},
                          true);

    io::pnetcdf::set_att(file,"averaging_window_size","NC_GLOBAL",m_time_avg_window_size[istream]);
  }

  io::pnetcdf::set_att(file,"start_date","NC_GLOBAL",m_case_t0.ymd());
  io::pnetcdf::set_att(file,"start_time","NC_GLOBAL",m_case_t0.tod());

  for (const auto& it1 : m_fields_stats[istream]) {
    for (const auto& it2 : it1.second) {
      const auto& stat  = it2.second;
      const auto& stat_layout = stat.layout();
      const auto& stat_dims = stat_layout.dims();
      const auto& stat_names = stat_layout.names();
      for (int axis = 0; axis < stat_layout.rank(); ++axis) {
        // TODO: this hard codes an E3SM impl detail (MPI decomposition over ncols)
        //       Add a setter method for host app to tell cldera which dim is decomposed
        //       over MPI.
        const bool decomposed = stat_names[axis]=="ncol";
        int dim = stat_dims[axis];
        io::pnetcdf::add_dim(file, stat_names[axis], dim, decomposed);
      }

      std::string io_dtype;
      if (stat.data_type()==DataType::RealType) {
        io_dtype = io::pnetcdf::get_io_dtype_name<Real>();
      } else if (stat.data_type()==DataType::IntType) {
        io_dtype = io::pnetcdf::get_io_dtype_name<int>();
      } else {
        EKAT_ERROR_MSG ("[ProfilingArchive::setup_output_file] Unsupported data type for IO.\n"
                        "  - stat name: " + stat.name() + "\n"
                        "  - data type: " + e2str(stat.data_type()) + "\n");
      }
      io::pnetcdf::add_var (file,
                            stat.name(),
                            io_dtype,
                            stat.layout().names(),
                            true);
    }
  }

  // List of fields (not stats) that we may need to write.
  // These are geometry-dep fields, like lat, lon, proc-rank,...
  std::list<std::string> non_stat_fields_to_write;

  // If requested, save geometry-related fields (if present)
  if (m_params.get<bool>("Save Geometry Fields",true)) {
    for (const auto& n : {"lat","lon","area"}) {
      if (has_field(n)) {
        // TODO: this hard codes an E3SM impl detail (MPI decomposition over ncols)
        //       Add a setter method for host app to tell cldera which dim is decomposed
        //       over MPI.
        io::pnetcdf::add_dim (file, "ncol", get_field(n).layout().size(),true);
        io::pnetcdf::add_var (file, n, io::pnetcdf::get_io_dtype_name<Real>(),{"ncol"},false);
        non_stat_fields_to_write.push_back(n);
      }
    }
  }

  // If any of the registered fields is distributed over columns,
  // we need to register the column MPI decomposition
  bool needs_decomp = false;
  if (file.dims.count("ncol")==1) {
    needs_decomp = true;

    // If we decide to split the output in N files (for size purposes), we need to
    // be careful, and not add the proc_rank field to the database more than once
    if (not has_field("proc_rank")) {
      // Also add the 'proc_rank' field, containing the owner of each col
      FieldLayout l ({file.dims.at("ncol")->len},{"ncol"});
      Field proc_rank("proc_rank",l,DataAccess::Copy,DataType::IntType);
      proc_rank.commit();
      add_field(proc_rank);
      Kokkos::deep_copy(proc_rank.view_nonconst<int>(),m_comm.rank());
    }
    io::pnetcdf::add_var (file, "proc_rank", io::pnetcdf::get_io_dtype_name<int>(),{"ncol"},false);

    io::pnetcdf::add_var (file, "col_gids", io::pnetcdf::get_io_dtype_name<int>(),{"ncol"},false);

    non_stat_fields_to_write.push_back("proc_rank");
    non_stat_fields_to_write.push_back("col_gids");
  }

  io::pnetcdf::enddef (file);
  if (needs_decomp) {
    auto f = get_field("col_gids");
    EKAT_REQUIRE_MSG (f.layout().rank()==1,
        "Error! Wrong layout for field 'col_gids'.\n"
        "  - expected layout: (ncol)\n"
        "  - actual layout  : (" + ekat::join(f.layout().names(),",") + ")\n");

    // Compute min gid. This is 1 in E3SM, but just in case...
    auto min_stat = StatFactory::instance().create("global_min",m_comm,ekat::ParameterList("global_min"));
    min_stat->set_field(f);
    min_stat->create_stat_field();
    auto min_gid = min_stat->compute(m_case_t0).data<int>()[0];

    int my_ngids = f.layout().size();
    std::vector<int> gids(my_ngids);
    for (int p=0,i=0; p<f.nparts(); ++p) {
      const auto pl = f.part_layout(p);
      const auto pdata = f.part_data<int>(p);
      for (int j=0; j<pl.size(); ++j,++i) {
        gids[i] = pdata[j] - min_gid;
      }
    }
    io::pnetcdf::add_decomp (file,"ncol",gids);
  }

  // Immediately write the non-time dep fields
  for (const auto& n : non_stat_fields_to_write) {
    const auto& f = get_field(n);
    // f might be partitioned, but our IO interface needs a single pointer to
    // the whole array. Therefore, if partitioned, we create a FieldIdentity
    // stat on the fly, compute it, and use the resulting 'stat' to pass the
    // data to the IO interface. Since this "stat" is computed only once, when
    // we setup the output file, the cost is tiny compared to the whole run.
    const int np = f.nparts();
    auto I = StatFactory::instance().create("identity",m_comm,ekat::ParameterList("identity"));
    I->set_field(f);
    I->create_stat_field();
    const auto f1p = np==1 ? f : I->compute(m_case_t0);
    if (f.data_type()==DataType::RealType) {
      io::pnetcdf::write_var (file,n,f1p.data<Real>());
    } else if (f.data_type()==DataType::IntType) {
      io::pnetcdf::write_var (file,n,f1p.data<int>());
    }
  }

  if (m_comm.am_i_root()) {
    printf("done!\n");
  }
  ts.stop_timer("profiling::setup_output_file");
}

void ProfilingArchive::
add_field (const Field& field)
{
  const auto& name = field.name();
  EKAT_REQUIRE_MSG (m_fields.find(name)==m_fields.end(),
      "[ProfilingArchive::add_field]\n"
      "  Error! Field '" + name + "' was already added.\n");

  m_fields_names.push_back(field.name());
  m_fields.emplace(name,field);
}

const Field& ProfilingArchive::
get_field (const std::string& name) const
{
  EKAT_REQUIRE_MSG (has_field(name),
      "[ProfilingArchive::get_field] Error! Field '" + name + "' not found.\n"
      "  List of current fields: " + ekat::join(m_fields_names,", "));
  return m_fields.at(name);
}

Field& ProfilingArchive::
get_field (const std::string& name)
{
  EKAT_REQUIRE_MSG (has_field(name),
      "[ProfilingArchive::get_field] Error! Field '" + name + "' not found.\n"
      "  List of current fields: " + ekat::join(m_fields_names,", "));
  return m_fields.at(name);
}

void ProfilingArchive::
update_stat (const std::string& fname, const std::string& stat_name,
             const Field& stat)
{
  EKAT_REQUIRE_MSG (has_field(fname),
      "[ProfilingArchive::get_field] Error! Field '" + fname + "' not found.\n"
      "  List of current fields: " + ekat::join(m_fields_names,", "));

  for (int i=0; i<m_num_streams; ++i) {
    auto& s = m_fields_stats[i][fname][stat_name];
    if (not s.committed()) {
      // It must be the first time we call update_stat for this stat.
      // Proceed to create the field
      s = stat.clone();
    }

    if (m_time_avg_curr_count[i]==0) {
      s.deep_copy(stat);
    } else {
      s.update(stat,1.0,1.0);
    }
  }
}

void ProfilingArchive::end_timestep (const TimeStamp& ts) {
  for (int i=0; i<m_num_streams; ++i) {
    ++m_time_avg_curr_count[i];
    if (m_time_avg_curr_count[i]==m_time_avg_window_size[i]) {
      // We completed averaging (or window size is 1, meaning no averaging), so
      //  1. set end of time avg window
      //  2. write vars
      //  3. reset beg of time avg window
      m_time_avg_end[i] = ts;
      write_stream(i);
      m_time_avg_beg[i] = ts;
      m_time_avg_curr_count[i] = 0;
    }
  }
}

void ProfilingArchive::write_stream (const int istream)
{
  auto& timings = timing::TimingSession::instance();
  timings.start_timer("profiling::write_stream");

  if (m_comm.am_i_root()) {
    printf(" [CLDERA] Flushing field stats to file ...\n");
  }

  auto& f = *m_output_files[istream];

  if (not f.enddef) {
    // We have not setup the output file yet
    setup_output_file(istream);
  }

  for (auto& it1 : m_fields_stats[istream]) {
    const auto& fname = it1.first;
    for (auto& it2: it1.second) {
      const auto& sname = it2.first;
            auto& stat  = it2.second;

      if (m_time_avg_window_size[istream]>1) {
        stat.scale (1.0/m_time_avg_window_size[istream]);
      }

      const auto& var_name = stat.name();

      if (stat.data_type()==DataType::RealType) {
        io::pnetcdf::write_var (f,var_name,stat.data<Real>());
      } else if (stat.data_type()==DataType::IntType) {
        io::pnetcdf::write_var (f,var_name,stat.data<int>());
      } else {
        EKAT_ERROR_MSG ("[ProfilingArchive::write_stream] Unsupported data type for IO.\n"
                        "  - stat name: " + stat.name() + "\n"
                        "  - data type: " + e2str(stat.data_type()) + "\n");
      }
    }
  }

  const double beg = m_time_avg_beg[istream] - m_case_t0;
  const double end = m_time_avg_end[istream] - m_case_t0;
  if (m_time_avg_window_size[istream]>1) {
    double bnds[2] = {beg,end};
    io::pnetcdf::write_var (f,"time_bounds",bnds);

    // We store beg as time, since that's when the avg window starts
    io::pnetcdf::update_time(f,beg);
  } else {
    io::pnetcdf::update_time(f,end);
  }

  if (m_comm.am_i_root()) {
    printf(" [CLDERA] Flushing field stats to file ... done!\n");
  }
  timings.stop_timer("profiling::write_stream");
}

void ProfilingArchive::commit_all_fields ()
{
  for (auto& it : m_fields) {
    it.second.commit();
  }
}

} // namespace cldera
