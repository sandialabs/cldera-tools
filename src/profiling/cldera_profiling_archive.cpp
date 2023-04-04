#include "cldera_profiling_archive.hpp"
#include "profiling/stats/cldera_field_stat_factory.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"
#include "timing/cldera_timing_session.hpp"

#include <ekat/io/ekat_yaml.hpp>
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <numeric>

namespace cldera {

ProfilingArchive::
ProfilingArchive(const ekat::Comm& comm,
                 const TimeStamp& start_date,
                 const ekat::ParameterList& params)
 : m_comm (comm)
 , m_params (params)
 , m_start_date (start_date)
{
  if (m_params.get<bool>("Enable Output",true)) {
    create_output_file();
  }
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
create_output_file ()
{
  const auto& file_name = m_params.get<std::string>("Filename","cldera_stats.nc");

  m_output_file = open_file (file_name, m_comm, io::pnetcdf::IOMode::Write);
  m_flush_freq = m_params.get("Flush Frequency",10);
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

  io::pnetcdf::set_att(*m_output_file,"start_date","NC_GLOBAL",m_start_date.ymd());
  io::pnetcdf::set_att(*m_output_file,"start_time","NC_GLOBAL",m_start_date.tod());

  for (const auto& it1 : m_fields_stats.front()) {
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
    auto min_stat = create_stat(ekat::ParameterList("global_min"),m_comm);
    auto min_gid = min_stat->compute(f).data<int>()[0];

    int my_ngids = f.layout().size();
    int ngids_scan;
    // Clock MPI ops
    track_mpi_scan(m_comm,&my_ngids,&ngids_scan,1,MPI_SUM,"profiling::setup_output_file");
    int my_start = ngids_scan - my_ngids;
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
    auto I = create_stat(ekat::ParameterList("identity"),m_comm);
    const auto f1p = np==1 ? f : I->compute(f);
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

  if (m_fields_stats.size()==0)
    EKAT_REQUIRE_MSG (m_time_stamps.size()==0,
        "Error! You had time stamps stored, but not stats.\n"
        "       You should have already gotten an error.\n"
        "       Please, contact developers.\n");

  if (m_fields_stats.size()==m_time_stamps.size())
    m_fields_stats.emplace_back();
  m_fields_stats.back()[fname][stat_name] = stat;
}

void ProfilingArchive::update_time (const TimeStamp& ts) {
  EKAT_REQUIRE_MSG (m_fields_stats.size()==(m_time_stamps.size()+1),
      "Error! It looks like you haven't stored any stat this time stamp.\n");
  m_time_stamps.push_back(ts);

  if (static_cast<int>(m_time_stamps.size())>m_flush_freq) {
    flush_to_file();
  }
}

void ProfilingArchive::flush_to_file ()
{
  if (m_output_file!=nullptr) {
    const int num_steps = m_time_stamps.size();
    if (num_steps == 0) {
      return;
    }

    auto& ts = timing::TimingSession::instance();
    ts.start_timer("profiling::flush_to_file");

    if (m_comm.am_i_root()) {
      printf(" [CLDERA] Flushing field stats to file ...\n");
    }

    if (not m_output_file->enddef) {
      // We have not setup the output file yet
      setup_output_file();
    }

    // Loop over number of time records
    auto ts_it = m_time_stamps.begin();
    auto stats_it = m_fields_stats.begin();
    for (int i=0; i<num_steps; ++i) {
      const auto& ts = *ts_it;
      const auto& stats = *stats_it;

      for (const auto& it1 : stats) {
        const auto& fname = it1.first;
        for (const auto& it2: it1.second) {
          const auto& sname = it2.first;
          const auto& stat  = it2.second;

          const auto var_name = fname + "_" + sname;

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
      io::pnetcdf::update_time(*m_output_file,ts-m_start_date);

      std::advance(ts_it,1);
      std::advance(stats_it,1);
    }

    m_time_stamps.clear();
    m_fields_stats.clear();

    if (m_comm.am_i_root()) {
      printf(" [CLDERA] Flushing field stats to file ... done!\n");
    }
    ts.stop_timer("profiling::flush_to_file");
  }
}

void ProfilingArchive::commit_all_fields ()
{
  for (auto& it : m_fields) {
    it.second.commit();
  }
}

} // namespace cldera
