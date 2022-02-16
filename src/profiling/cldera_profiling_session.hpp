#ifndef CLDERA_PROFILING_SESSION_HPP
#define CLDERA_PROFILING_SESSION_HPP

#include <ekat/mpi/ekat_comm.hpp>
#include <ekat/std_meta/ekat_std_any.hpp>

#include <string>
#include <map>

namespace cldera {

/*
 * Singleton class holding all the profiling data
 *
 * The class is storing a map string->ekat::any.
 * There's zero knowledge of what data is. The goal
 * of this class is just to provide a centralized,
 * persistent storage for data, so that the profiling
 * library may be called from an external app,
 * without the need for the app to keep track of
 * all the cldera data structures.
 *
 * See ekat_std_any.hpp for how to create an ekat::any
 * object, and for how to retrieved the internal data.
 */

class ProfilingSession {
public:
  static ProfilingSession& instance () {
    static ProfilingSession s;
    return s;
  };

  void init (const ekat::Comm& comm);
  void clean_up ();

  const ekat::Comm& get_comm () const {
    return m_comm;
  }

  void set_data (const std::string& name, ekat::any& data);
  const ekat::any& get_data (const std::string& name) const;

private:
  ProfilingSession () = default;

  std::map<std::string,ekat::any> m_data;

  ekat::Comm  m_comm;
  bool m_inited = false;
};

} // namespace cldera

#endif // CLDERA_PROFILING_SESSION_HPP
