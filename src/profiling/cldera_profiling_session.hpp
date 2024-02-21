#ifndef CLDERA_PROFILING_SESSION_HPP
#define CLDERA_PROFILING_SESSION_HPP

#include "cldera_profiling_context.hpp"

#include <map>
#include <string>

namespace cldera
{

/*
 * A ProfilingContext is a structure meant to hold multiple
 * ProfilingSessions objects. The user needs to change the
 * name of the session they want to access.
 */
struct ProfilingSession {
  static ProfilingSession& instance () {
    static ProfilingSession s;
    return s;
  }

  void init_session (const ekat::Comm& comm);

  ProfilingContext& add_context (const std::string& name);

  ProfilingContext& get_curr_context ();

  void switch_context (const std::string& name);

  void clean_up_current_context ();

  const std::string& curr_context_name () const { return m_curr_context_name; }

private:

  std::map<std::string,ProfilingContext> m_contexts;

  ProfilingSession () = default;

  std::string   m_curr_context_name;
};

} // namespace cldera

#endif // CLDERA_PROFILING_SESSION_HPP
