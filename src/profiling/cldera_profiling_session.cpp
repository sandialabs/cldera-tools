#include "cldera_profiling_session.hpp"
#include "cldera_version.h"

#include <ekat/ekat_session.hpp>

namespace cldera
{

void ProfilingSession::
init_session (const ekat::Comm& comm)
{
  static bool ekat_inited = false;
  if (not ekat_inited) {
    ekat::initialize_ekat_session (/* argc = */ 0,
                                   /* argv = */ nullptr,
                                   /* print config = */ comm.am_i_root());
    ekat_inited = true;
    if (comm.am_i_root()) {
      std::cout << "  -> CLDERA version: " << CLDERA_VERSION << "\n"
                   "  -> CLDERA git sha: " << CLDERA_GIT_SHA << "\n";
    }
  }
}

ProfilingContext& ProfilingSession::
add_context (const std::string& name)
{
  EKAT_REQUIRE_MSG (m_contexts.count(name)==0,
      "Error! ProfilingSession was already inited.\n"
      "  name: " + name + "\n");

  m_curr_context_name = name;
  auto& c = m_contexts.emplace(name,name).first->second;
  return c;
}

ProfilingContext&
ProfilingSession::get_curr_context ()
{
  EKAT_REQUIRE_MSG (m_contexts.count(m_curr_context_name)==1,
      "Error! No ProfilingContext currently active.\n");
  return m_contexts.at(m_curr_context_name);
}

void ProfilingSession::
switch_context (const std::string& name)
{
  EKAT_REQUIRE_MSG (m_contexts.count(name)==1,
      "Error! ProfilingContext was never create.\n"
      "  name: " + name + "\n");
  m_curr_context_name = name;
}

void ProfilingSession::
clean_up_current_context ()
{
  // Clean up and erase current context
  get_curr_context().clean_up();
  m_contexts.erase(m_curr_context_name);

  // No "current" context anymore
  m_curr_context_name = "";

  // If this was the last context, shut down ekat
  if (m_contexts.size()==0) {
    ekat::finalize_ekat_session ();
  }
}

} // namespace cldera
