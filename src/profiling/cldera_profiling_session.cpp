#include "profiling/cldera_profiling_session.hpp"

#include <ekat/ekat_assert.hpp>

namespace cldera {

void ProfilingSession::
init (const ekat::Comm& comm)
{
  EKAT_REQUIRE_MSG (not m_inited, "Error! ProfilingSession was already inited.\n");
  m_comm = comm;
  m_inited = true;
}

void ProfilingSession::
clean_up ()
{
  EKAT_REQUIRE_MSG (m_inited, "Error! ProfilingSession was not yet inited.\n");
  m_comm.reset_mpi_comm(MPI_COMM_NULL);
  m_inited = false;
}

void ProfilingSession::
set_data (const std::string& name, ekat::any& data)
{
  EKAT_REQUIRE_MSG (m_data.find(name)==m_data.end(),
      "[cldera::ProfilingSession] Error! Data '" + name + "' was already set.\n");

  m_data[name] = data;
}

const ekat::any& ProfilingSession::
get_data (const std::string& name) const
{
  auto it = m_data.find(name);
  EKAT_REQUIRE_MSG (it!=m_data.end(),
      "[cldera::ProfilingSession] Error! Data '" + name + "' not found.\n");

  return it->second;
}

} // namespace cldera
