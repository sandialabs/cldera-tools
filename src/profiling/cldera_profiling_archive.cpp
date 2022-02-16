#include "cldera_profiling_archive.hpp"

#include <ekat/ekat_assert.hpp>

namespace cldera {

void ProfilingArchive::
clean_up () {
  m_fields.clear();
  m_field_stats.clear();
}

void ProfilingArchive::
add_field (const std::string& name,
           const Real* const data,
           const std::vector<int>& dims)
{
  add_field (name,Field{data,dims});
}

void ProfilingArchive::
add_field (const std::string& name,
           const Field& field)
{
  EKAT_REQUIRE_MSG (m_fields.find(name)==m_fields.end(),
      "[ProfilingArchive::add_field] Error! Field '" + name + "' was already added.\n");

  EKAT_REQUIRE_MSG (field.data!=nullptr,
      "[ProfilingArchive::add_field] Error! Input Field '" + name + "' stores a nullptr.\n");

  m_fields.emplace(name,field);
  m_field_stats.emplace(name,stats_history_t{});

  // Ensure a (empty) History is present for all stat types.
  stats_history_t stats;
  stats.emplace(StatType::Avg, History{});
  stats.emplace(StatType::Max, History{});
  stats.emplace(StatType::Min, History{});
  stats.emplace(StatType::Sum, History{});
}

const Field& ProfilingArchive::
get_field (const std::string& name) const
{
  auto it = m_fields.find(name);
  EKAT_REQUIRE_MSG (it!=m_fields.end(),
      "[ProfilingArchive::get_field] Error! Field '" + name + "' not found.\n");

  return it->second;
}

const ProfilingArchive::stats_history_t&
ProfilingArchive::
get_all_stats_history (const std::string& name) const
{
  auto it = m_field_stats.find(name);
  EKAT_REQUIRE_MSG (it!=m_field_stats.end(),
      "[ProfilingArchive::get_field] Error! Stats for field '" + name + "' not found.\n");

  return it->second;
}

const History&
ProfilingArchive::
get_stat_history (const std::string& name, const StatType stat) const
{
  return get_all_stats_history(name).at(stat);
}

void ProfilingArchive::
store_stat (const std::string& name,
            const StatType     stat_type,
            const Real         time,
            const Real         value)
{
  auto it = m_field_stats.find(name);
  EKAT_REQUIRE_MSG (it!=m_field_stats.end(),
      "[ProfilingArchive::add_field] Error! Field '" + name + "' was not registered.\n");
  auto& history = it->second[stat_type];

  EKAT_REQUIRE_MSG (history.times.size()==0 || time>history.times.back(),
      "[ProfilingArchive::add_field] Error! Time stamp for stat is not increasing.\n"
      "   - Field name : " + name + ".\n"
      "   - Last time  : " + std::to_string(history.times.back()) + ".\n"
      "   - Input time : " + std::to_string(value) + ".\n");

  history.times.push_back(time);
  history.values.push_back(value);
}

} // namespace cldera
