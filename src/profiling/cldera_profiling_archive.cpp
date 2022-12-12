#include "cldera_profiling_archive.hpp"
#include "profiling/stats/cldera_field_stat_factory.hpp"

#include <ekat/io/ekat_yaml.hpp>
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera {

ProfilingArchive::
ProfilingArchive(const ekat::Comm& comm,
                 const TimeStamp& start_date,
                 const ekat::ParameterList& params)
 : m_comm (comm)
 , m_start_date (start_date)
{
  create_output_file(params);
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
create_output_file (const ekat::ParameterList& params)
{
  if (params.isParameter("Filename")) {
    const auto& file_name = params.get<std::string>("Filename");
    m_output_file = open_file (file_name, m_comm, io::pnetcdf::IOMode::Write);
  } else {
    const std::string file_name = "cldera_stats.nc";
    m_output_file = open_file (file_name, m_comm, io::pnetcdf::IOMode::Write);
  }
  if (params.isParameter("Flush Frequency")) {
    m_flush_freq = params.get<int>("Flush Frequency");
  } else {
    m_flush_freq = 10;
  }
}

void ProfilingArchive::
setup_output_file ()
{
  if (m_comm.am_i_root()) {
    printf(" [CLDERA] setting up output file ...");
  }

  // Add time dimension
  io::pnetcdf::add_time (*m_output_file,"double");

  io::pnetcdf::set_att(*m_output_file,"start_date","NC_GLOBAL",m_start_date.ymd());
  io::pnetcdf::set_att(*m_output_file,"start_time","NC_GLOBAL",m_start_date.tod());

  for (const auto& it1 : m_fields_stats.front()) {
    const auto& fname = it1.first;
    const auto& f = get_field(fname);
    for (const auto& it2 : it1.second) {
      const auto& sname = it2.first;
      const auto& stat  = it2.second;
      const auto& stat_layout = stat.layout();
      const auto& stat_dims = stat_layout.dims();
      const auto& stat_names = stat_layout.names();
      for (int axis = 0; axis < stat_layout.rank(); ++axis)
        add_dim(*m_output_file, stat_names[axis], stat_dims[axis]);
      add_var (*m_output_file,
               stat.name(),
               io::pnetcdf::get_io_dtype_name<Real>(),
               stat.layout().names(),
               true);
    }
  }
  io::pnetcdf::enddef (*m_output_file);

  if (m_comm.am_i_root()) {
    printf("done!\n");
  }
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

    if (m_comm.am_i_root()) {
      printf(" [CLDERA] Flushing field stats to file ...\n");
    }

    const int num_steps = m_time_stamps.size();
    if (num_steps == 0) {
      return;
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

          write_var (*m_output_file,var_name,stat.data<Real>());
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
  }
}

void ProfilingArchive::commit_all_fields ()
{
  for (auto& it : m_fields) {
    it.second.commit();
  }
}

} // namespace cldera
