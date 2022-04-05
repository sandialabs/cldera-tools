#include "cldera_profiling_archive.hpp"

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

  m_fields.emplace(name,field);
}

const Field& ProfilingArchive::
get_field (const std::string& name) const
{
  EKAT_REQUIRE_MSG (has_field(name),
      "[ProfilingArchive::get_field] Error! Field '" + name + "' not found.\n");
  return m_fields.at(name);
}

Field& ProfilingArchive::
get_field (const std::string& name)
{
  EKAT_REQUIRE_MSG (has_field(name),
      "[ProfilingArchive::get_field] Error! Field '" + name + "' not found.\n");
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

void ProfilingArchive::commit_all_fields ()
{
  for (auto& it : m_fields) {
    it.second.commit();
  }
}
} // namespace cldera
