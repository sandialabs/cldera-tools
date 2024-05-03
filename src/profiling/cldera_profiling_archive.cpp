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
  if (m_params.get<bool>("Enable Output",true)) {
    const auto& prefix = m_params.get<std::string>("filename_prefix","cldera_stats");
    std::string filename;
    io::pnetcdf::IOMode mode = io::pnetcdf::IOMode::Invalid;

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
    m_output_file = open_file (filename+".nc", m_comm, mode);
  }
  m_flush_freq = m_params.get("Flush Frequency",10);

  // Init the beg/end vector for time averaging
  m_time_stamps_beg.resize(m_flush_freq,run_t0);
  m_time_stamps_end.resize(m_flush_freq);

  m_time_avg_window_size = m_params.get("Time Averaging Window Size",1);
}

ProfilingArchive::
~ProfilingArchive()
{
  // If there's any data left to flush, do it.
  flush_to_file ();

  if (m_output_file) {
    io::pnetcdf::close_file(*m_output_file);
  }
}

void ProfilingArchive::
setup_output_file ()
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("profiling::setup_output_file");

  if (m_comm.am_i_root()) {
    printf(" [CLDERA] setting up output file ...");
  }

  // Add time dimension
  io::pnetcdf::add_time (*m_output_file,"double");

  if (m_time_avg_window_size>1) {
    io::pnetcdf::add_dim(*m_output_file, "dim2", 2);
    // We save the start/end of each averaging window
    io::pnetcdf::add_var (*m_output_file,
                           "time_bounds",
                           io::pnetcdf::get_io_dtype_name<Real>(),
                           {"dim2"},
                           true);

    io::pnetcdf::set_att(*m_output_file,"averaging_window_size","NC_GLOBAL",m_time_avg_window_size);
  }

  io::pnetcdf::set_att(*m_output_file,"start_date","NC_GLOBAL",m_case_t0.ymd());
  io::pnetcdf::set_att(*m_output_file,"start_time","NC_GLOBAL",m_case_t0.tod());

  for (const auto& it1 : m_fields_stats) {
    for (const auto& it2 : it1.second) {
      const auto& stat  = it2.second[0];
      const auto& stat_layout = stat.layout();
      const auto& stat_dims = stat_layout.dims();
      const auto& stat_names = stat_layout.names();
      for (int axis = 0; axis < stat_layout.rank(); ++axis) {
        // TODO: this hard codes an E3SM impl detail (MPI decomposition over ncols)
        //       Add a setter method for host app to tell cldera which dim is decomposed
        //       over MPI.
        const bool decomposed = stat_names[axis]=="ncol";
        int dim = stat_dims[axis];
        io::pnetcdf::add_dim(*m_output_file, stat_names[axis], dim, decomposed);
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
      io::pnetcdf::add_var (*m_output_file,
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
        io::pnetcdf::add_dim (*m_output_file, "ncol", get_field(n).layout().size(),true);
        io::pnetcdf::add_var (*m_output_file, n, io::pnetcdf::get_io_dtype_name<Real>(),{"ncol"},false);
        non_stat_fields_to_write.push_back(n);
      }
    }
  }

  // If any of the registered fields is distributed over columns,
  // we need to register the column MPI decomposition
  bool needs_decomp = false;
  if (m_output_file->dims.count("ncol")==1) {
    needs_decomp = true;

    // If we decide to split the output in N files (for size purposes), we need to
    // be careful, and not add the proc_rank field to the database more than once
    if (not has_field("proc_rank")) {
      // Also add the 'proc_rank' field, containing the owner of each col
      FieldLayout l ({m_output_file->dims.at("ncol")->len},{"ncol"});
      Field proc_rank("proc_rank",l,DataAccess::Copy,DataType::IntType);
      proc_rank.commit();
      add_field(proc_rank);
      Kokkos::deep_copy(proc_rank.view_nonconst<int>(),m_comm.rank());
    }
    io::pnetcdf::add_var (*m_output_file, "proc_rank", io::pnetcdf::get_io_dtype_name<int>(),{"ncol"},false);

    io::pnetcdf::add_var (*m_output_file, "col_gids", io::pnetcdf::get_io_dtype_name<int>(),{"ncol"},false);

    non_stat_fields_to_write.push_back("proc_rank");
    non_stat_fields_to_write.push_back("col_gids");
  }

  io::pnetcdf::enddef (*m_output_file);
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
    io::pnetcdf::add_decomp (*m_output_file,"ncol",gids);
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
      io::pnetcdf::write_var (*m_output_file,n,f1p.data<Real>());
    } else if (f.data_type()==DataType::IntType) {
      io::pnetcdf::write_var (*m_output_file,n,f1p.data<int>());
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
append_stat (const std::string& fname, const std::string& stat_name,
             const Field& stat)
{
  EKAT_REQUIRE_MSG (has_field(fname),
      "[ProfilingArchive::get_field] Error! Field '" + fname + "' not found.\n"
      "  List of current fields: " + ekat::join(m_fields_names,", "));

  auto update_stat = [&](Field& curr_f, const Field& new_f) {
    if (m_time_avg_curr_count==0) {
      curr_f = new_f.clone();
    } else {
      curr_f.update(new_f,1,1);
    }
  };
  auto& vec = m_fields_stats[fname][stat_name];

  // Only done the first time
  vec.reserve(m_flush_freq);

  vec.resize(m_curr_time_slot);
  update_stat(vec.back(),stat);
}

void ProfilingArchive::end_timestep (const TimeStamp& ts) {
  ++m_time_avg_curr_count;
  if (m_time_avg_curr_count==m_time_avg_window_size) {
    // We completed averaging (or window size is 1, meaning no averagin),
    // so store end timestamp for this window, and set beg for next window
    m_time_stamps_end[m_curr_time_slot] = ts;
    ++m_curr_time_slot;
    if (m_curr_time_slot==m_flush_freq) {
      flush_to_file();
    }
    m_time_stamps_beg[m_curr_time_slot] = ts;
    m_time_avg_curr_count = 0;
  }
}

void ProfilingArchive::flush_to_file ()
{
  // Note: if m_curr_time_slot=0, then we are being called from the destructor,
  // and no additional stat has been computed since the last flush.
  if (m_output_file!=nullptr and m_curr_time_slot>0) {

    auto& timings = timing::TimingSession::instance();
    timings.start_timer("profiling::flush_to_file");

    if (m_comm.am_i_root()) {
      printf(" [CLDERA] Flushing field stats to file ...\n");
    }

    if (not m_output_file->enddef) {
      // We have not setup the output file yet
      setup_output_file();
    }

    // Loop over number of time records
    for (int i=0; i<m_curr_time_slot; ++i) {
      const double beg = m_time_stamps_beg[i] - m_case_t0;
      const double end = m_time_stamps_end[i] - m_case_t0;

      for (const auto& it1 : m_fields_stats) {
        const auto& fname = it1.first;
        for (const auto& it2: it1.second) {
          const auto& sname = it2.first;
          const auto& stat  = it2.second[i];

          const auto var_name = stat.name();

          if (stat.data_type()==DataType::RealType) {
            io::pnetcdf::write_var (*m_output_file,var_name,stat.data<Real>());
          } else if (stat.data_type()==DataType::IntType) {
            io::pnetcdf::write_var (*m_output_file,var_name,stat.data<int>());
          } else {
            EKAT_ERROR_MSG ("[ProfilingArchive::flush_to_file] Unsupported data type for IO.\n"
                            "  - stat name: " + stat.name() + "\n"
                            "  - data type: " + e2str(stat.data_type()) + "\n");
          }
        }
      }
      if (m_time_avg_window_size>1) {
        double bnds[2] = {beg,end};
        io::pnetcdf::write_var (*m_output_file,"time_bounds",bnds);

        // We store beg as time, since that's when the avg window starts
        io::pnetcdf::update_time(*m_output_file,beg);
      } else {
        io::pnetcdf::update_time(*m_output_file,end);
      }

    }

    if (m_comm.am_i_root()) {
      printf(" [CLDERA] Flushing field stats to file ... done!\n");
    }
    timings.stop_timer("profiling::flush_to_file");
  }

  // Whether we do have an output file or not, we must reset the time slot to 0
  m_curr_time_slot = 0;
}

void ProfilingArchive::commit_all_fields ()
{
  for (auto& it : m_fields) {
    it.second.commit();
  }
}

} // namespace cldera
