#ifndef CLDERA_PROFILING_ARCHIVE_HPP
#define CLDERA_PROFILING_ARCHIVE_HPP

#include "cldera_profiling_types.hpp"

#include <map>

namespace cldera {

/*
 * Storage for data to profile
 *
 * The storage includes both information on the arrays
 * to track/profile, as well as a history of their
 * computed statistics.
 *
 * This class is not computing anything. It is simply
 * acting as a repository for all the fields and their stats.
 */

class ProfilingArchive
{
public:
  using stats_history_t = std::map<StatType,History>;

  ~ProfilingArchive ();

  void clean_up ();

  // Setters
  void add_field (const std::string& name,
                  const Real* const data,
                  const std::vector<int>& dims);

  void add_field (const std::string& name,
                  const Field& field);

  void store_stat (const std::string& name,
                   const StatType     stat_type,
                   const Real         time,
                   const Real         value);

  // Getters
  const Field& get_field (const std::string& name) const;
  stats_history_t& get_all_stats_history (const std::string& name);
  History& get_stat_history (const std::string& name, const StatType stat);

private:

  std::map<std::string,Field>             m_fields;
  std::map<std::string,stats_history_t>   m_field_stats;
};

} // namespace cldera

#endif // CLDERA_PROFILING_ARCHIVE_HPP
