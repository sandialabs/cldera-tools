#ifndef CLDERA_PROFILING_CONTEXT_HPP
#define CLDERA_PROFILING_CONTEXT_HPP

#include "profiling/cldera_time_stamp.hpp"

#include "timing/cldera_timing_session.hpp"

#include <ekat/mpi/ekat_comm.hpp>
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/std_meta/ekat_std_any.hpp>

#include <string>
#include <map>

namespace cldera {

/*
 * Class holding data for profiling
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

class ProfilingContext {
public:
  ProfilingContext (const std::string& name)
    : m_name(name) {}

  void init (const ekat::Comm& comm,
             const ekat::ParameterList& params);

  // Note: you may think to do this stuff in the destructor. However,
  //       there is a chance this call may throw, and you CAN'T throw
  //       from inside a destructor
  void clean_up ();

  const ekat::Comm& get_comm () const { return m_comm; }

  const ekat::ParameterList& get_params () const { return m_params; }
        ekat::ParameterList& get_params ()       { return m_params; }

  bool has_data (const std::string& name) const {
    return m_data.find(name)!=m_data.end();
  }

  // Construct the any object on the fly
  template<typename T, typename... Args>
  T& create (const std::string& name, const Args&... args);

  template<typename T>
  T& get (const std::string& name) const;

  bool inited () const { return m_inited; }

  const std::string& name () const { return m_name; }

  timing::TimingSession& timing () { return m_timing; }

private:

  timing::TimingSession m_timing;

  std::map<std::string,ekat::any> m_data;

  ekat::Comm  m_comm;
  ekat::ParameterList m_params;
  bool m_inited = false;

  std::string m_name;
};

// ============================= IMPLEMENTATION ============================= //

// Construct the any object on the fly
template<typename T, typename... Args>
T& ProfilingContext::
create (const std::string& name, const Args&... args)
{
  EKAT_REQUIRE_MSG (not has_data(name),
      "[cldera::ProfilingContext] Error! Data '" + name + "' was already created.\n");

  m_data[name].reset<T>(args...);

  return get<T>(name);
}

template<typename T>
T& ProfilingContext::get (const std::string& name) const
{
  EKAT_REQUIRE_MSG (has_data(name),
      "[cldera::ProfilingContext] Error! Data '" + name + "' not found.\n");

  auto data = m_data.at(name);

  return ekat::any_cast<T>(data);
}

} // namespace cldera

#endif // CLDERA_PROFILING_CONTEXT_HPP
