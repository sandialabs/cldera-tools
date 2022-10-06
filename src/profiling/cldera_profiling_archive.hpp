#ifndef CLDERA_PROFILING_ARCHIVE_HPP
#define CLDERA_PROFILING_ARCHIVE_HPP

#include "profiling/stats/cldera_field_stat.hpp"
#include "cldera_profiling_types.hpp"
#include "cldera_time_stamp.hpp"
#include "cldera_field.hpp"

#include <io/cldera_pnetcdf.hpp>

#include <ekat/ekat_parameter_list.hpp>

#include <map>
#include <list>

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
  template<typename T>
  using strmap_t = std::map<std::string,T>;

  using stat_ptr_t = std::shared_ptr<FieldStat>;

  ProfilingArchive (const ekat::Comm& comm,
                    const TimeStamp& start_date,
                    const ekat::ParameterList& params);

  // Flushes data to file, if needed
  ~ProfilingArchive ();

  // Fields
  void add_field (const Field& field);

  void commit_all_fields ();

  bool has_field (const std::string& name) const {
    return m_fields.find(name)!=m_fields.end();
  }
  const Field& get_field (const std::string& name) const;
        Field& get_field (const std::string& name);

  // Stats
  void append_stat (const std::string& fname, const std::string& stat_name,
                    const Field& stat);

  void update_time (const TimeStamp& ts);
private:
  void create_output_file (const ekat::ParameterList& params);
  void setup_output_file ();

  void flush_to_file ();

  ekat::Comm                              m_comm;

  std::list<std::string>                  m_fields_names;

  strmap_t<Field>                         m_fields;
  std::list<strmap_t<strmap_t<Field>>>    m_fields_stats;

  TimeStamp                               m_start_date;
  std::list<TimeStamp>                    m_time_stamps;

  std::shared_ptr<io::pnetcdf::NCFile>    m_output_file;

  // TODO: make this a runtime option
  int m_flush_freq = 10;
};

} // namespace cldera

#endif // CLDERA_PROFILING_ARCHIVE_HPP
