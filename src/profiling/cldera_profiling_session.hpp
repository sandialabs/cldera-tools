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

  bool has_data (const std::string& name) const {
    return m_data.find(name)!=m_data.end();
  }

  // Construct the any object on the fly
  template<typename T, typename... Args>
  T& create (const std::string& name, const Args&... args);

  // If not already there, create it
  template<typename T, typename... Args>
  T& create_or_get(const std::string& name, const Args&... args);

  template<typename T>
  T& get (const std::string& name) const;

  void remove (const std::string& name);

  bool inited () const { return m_inited; }
private:
  ProfilingSession () = default;

  std::map<std::string,ekat::any> m_data;

  ekat::Comm  m_comm;
  bool m_inited = false;
};

// ============================= IMPLEMENTATION ============================= //

// Construct the any object on the fly
template<typename T, typename... Args>
T& ProfilingSession::
create (const std::string& name, const Args&... args)
{
  EKAT_REQUIRE_MSG (not has_data(name),
      "[cldera::ProfilingSession] Error! Data '" + name + "' was already created.\n");

  m_data[name].reset<T>(args...);

  return get<T>(name);
}

template<typename T, typename... Args>
T& ProfilingSession::
create_or_get (const std::string& name, const Args&... args)
{
  if (not has_data(name)) {
    return create<T>(name,args...);
  }

  return get<T>(name);
}

template<typename T>
T& ProfilingSession::get (const std::string& name) const
{
  EKAT_REQUIRE_MSG (has_data(name),
      "[cldera::ProfilingSession] Error! Data '" + name + "' not found.\n");

  auto data = m_data.at(name);

  return ekat::any_cast<T>(data);
}

} // namespace cldera

#endif // CLDERA_PROFILING_SESSION_HPP
