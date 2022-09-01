#include "cldera_profiling_archive.hpp"
#include "stats/cldera_stats_utils.hpp"

#include <ekat/io/ekat_yaml.hpp>
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera {

ProfilingArchive::~ProfilingArchive()
{
  clean_up();
}

void ProfilingArchive::
clean_up () {
  m_fields.clear();
  m_field_stats.clear();
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

History&
ProfilingArchive::
get_stat_history (const std::string& name, const StatType stat)
{
  EKAT_REQUIRE_MSG (has_field(name),
      "[ProfilingArchive::get_field] Error! Field '" + name + "' not registered.\n");

  return m_field_stats[name][stat];
}

void ProfilingArchive::dump_stats_to_yaml (const std::string& file_name) const
{
  ekat::ParameterList stats("Statistics");

  // Loop over fields
  for (const auto& it1 : m_field_stats) {
    const auto& fname = it1.first;

    // Loop over all stats for this field
    for (const auto& it2 : it1.second) {
      auto stat = it2.first;
      const auto& hist = it2.second;

      // If hist is empty, skip it
      if (hist.size()>0) {
        auto& stats_of_f = stats.sublist(fname);
        auto& this_stat = stats_of_f.sublist(e2str(stat));
        this_stat.set("Values",hist.values());
        auto& times = this_stat.get<std::vector<std::string>>("Timestamps",{});

        // Write timestamps as YYYYMMDD.TOD
        for (const auto& t : hist.times()) {
          times.push_back(std::to_string(t.ymd) + "." + std::to_string(t.tod));
        }
      }
    }
  }

  ekat::write_yaml_file(file_name,stats);
}

void ProfilingArchive::commit_all_fields ()
{
  for (auto& it : m_fields) {
    it.second.commit();
  }
}
} // namespace cldera
