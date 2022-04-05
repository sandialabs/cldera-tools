#ifndef CLDERA_PROFILING_ARCHIVE_HPP
#define CLDERA_PROFILING_ARCHIVE_HPP

#include "cldera_profiling_types.hpp"
#include "cldera_field.hpp"

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

  // Fields
  void add_field (const Field& field);

  void commit_all_fields ();

  bool has_field (const std::string& name) const {
    return m_fields.find(name)!=m_fields.end();
  }
  const Field& get_field (const std::string& name) const;
        Field& get_field (const std::string& name);

  // Stats
  History& get_stat_history (const std::string& name, const StatType stat);
private:

  std::map<std::string,Field>             m_fields;
  std::map<std::string,stats_history_t>   m_field_stats;
};

} // namespace cldera

#endif // CLDERA_PROFILING_ARCHIVE_HPP
