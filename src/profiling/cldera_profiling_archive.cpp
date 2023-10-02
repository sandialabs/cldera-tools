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
    m_output_file = open_file (filename, m_comm, mode);
  }
  m_flush_freq = m_params.get("Flush Frequency",10);
  m_time_stamps.resize(m_flush_freq);
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
        add_dim(*m_output_file, stat_names[axis], dim, decomposed);
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
      add_var (*m_output_file,
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

  auto& vec = m_fields_stats[fname][stat_name];
  if (vec.size()==0) {
    // It's the first time we append to this stat. Create the vector
    vec.resize(m_flush_freq);
    for (auto& f : vec)
      f = stat.clone();
  } else {
    EKAT_REQUIRE_MSG (m_curr_time_slot<vec.size(),
        "Error! We ran out of slots for this stat. Looks like you missed a call to flush_to_file?\n"
        " - Field name: " + fname + "\n"
        " - Stat name : " + stat_name + "\n");
    vec[m_curr_time_slot].deep_copy(stat);
  }
}

void ProfilingArchive::update_time (const TimeStamp& ts) {
  m_time_stamps[m_curr_time_slot] = ts;

  ++m_curr_time_slot;
  if (m_curr_time_slot==m_flush_freq) {
    flush_to_file();
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
      const auto& ts = m_time_stamps[i];

      for (const auto& it1 : m_fields_stats) {
        const auto& fname = it1.first;
        for (const auto& it2: it1.second) {
          const auto& sname = it2.first;
          const auto& stat  = it2.second[i];

          const auto var_name = stat.name();

          if (stat.data_type()==DataType::RealType) {
            write_var (*m_output_file,var_name,stat.data<Real>());
          } else if (stat.data_type()==DataType::IntType) {
            write_var (*m_output_file,var_name,stat.data<int>());
          } else {
            EKAT_ERROR_MSG ("[ProfilingArchive::flush_to_file] Unsupported data type for IO.\n"
                            "  - stat name: " + stat.name() + "\n"
                            "  - data type: " + e2str(stat.data_type()) + "\n");
          }
        }
      }
      io::pnetcdf::update_time(*m_output_file,ts-m_case_t0);
    }

    m_curr_time_slot = 0;

    if (m_comm.am_i_root()) {
      printf(" [CLDERA] Flushing field stats to file ... done!\n");
    }
    timings.stop_timer("profiling::flush_to_file");
  }
}

void ProfilingArchive::commit_all_fields ()
{
  for (auto& it : m_fields) {
    it.second.commit();
  }
}

} // namespace cldera
