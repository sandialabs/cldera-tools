#include "profiling/cldera_profiling_session.hpp"

#include <ekat/ekat_assert.hpp>

namespace cldera {

void ProfilingSession::
init (const ekat::Comm& comm)
{
  EKAT_REQUIRE_MSG (not inited(), "Error! ProfilingSession was already inited.\n");
  m_comm = comm;
  m_inited = true;
}

void ProfilingSession::
clean_up ()
{
  EKAT_REQUIRE_MSG (inited(), "Error! ProfilingSession was not yet inited.\n");

  m_comm.reset_mpi_comm(MPI_COMM_SELF);
  m_data.clear();

  m_inited = false;
}

void ProfilingSession::
remove (const std::string& name)
{
  EKAT_REQUIRE_MSG (has_data(name),
      "[cldera::ProfilingSession] Error! Data '" + name + "' not stored.\n");

  m_data.erase(name);
}

} // namespace cldera
