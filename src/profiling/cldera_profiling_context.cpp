#include "profiling/cldera_profiling_context.hpp"
#include "profiling/cldera_profiling_archive.hpp"
#include "profiling/cldera_pathway_factory.hpp"
#include "cldera_pathway.hpp"

#include <ekat/ekat_assert.hpp>

#include <fstream>

namespace cldera {

void ProfilingContext::
init (const ekat::Comm& comm,
      const ekat::ParameterList& params)
{
  EKAT_REQUIRE_MSG (not inited(), "Error! ProfilingContext was already inited.\n");
  m_comm = comm;
  m_params = params;
  m_inited = true;

  // If no filename is provided, there's no point in doing timing
  m_timing.toggle_session(m_params.isParameter("Timing Filename"));
}

void ProfilingContext::
clean_up ()
{
  EKAT_REQUIRE_MSG (inited(), "Error! ProfilingContext was not yet inited.\n");
  m_timing.start_timer(m_name + "::clean_up");

  m_data.clear();

  m_timing.stop_timer(m_name + "::clean_up");

  if (m_timing.is_active()) {
    const auto& timing_fname = m_params.get<std::string>("Timing Filename");
    std::ofstream timing_file;
    std::stringstream blackhole;
    if (m_comm.am_i_root()) {
      timing_file.open(timing_fname);
    }
    std::ostream& ofile = timing_file;
    std::ostream& onull = blackhole;

    std::ostream& out = m_comm.am_i_root() ? ofile : onull;
    m_timing.dump(out,m_comm);
  }
  m_timing.clean_up();

  m_comm.reset_mpi_comm(MPI_COMM_SELF);
  m_inited = false;
}

} // namespace cldera
