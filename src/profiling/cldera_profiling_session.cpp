#include "profiling/cldera_profiling_session.hpp"

#include <ekat/ekat_assert.hpp>

namespace cldera {

void ProfilingSession::
init (const ekat::Comm& comm) {
  EKAT_REQUIRE_MSG (not m_inited, "Error! ProfilingSession was already inited.\n");
  m_comm = comm;
  m_inited = true;
}

void ProfilingSession::
clean_up () {
  EKAT_REQUIRE_MSG (m_inited, "Error! ProfilingSession was not yet inited.\n");
  m_comm.reset_mpi_comm(MPI_COMM_NULL);
  m_inited = false;
}

} // namespace cldera
