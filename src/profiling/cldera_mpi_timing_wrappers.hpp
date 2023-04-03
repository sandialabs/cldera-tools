#ifndef CLDERA_MPI_TIMING_WRAPPERS_HPP
#define CLDERA_MPI_TIMING_WRAPPERS_HPP

#include "timing/cldera_timing_session.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

template<typename T>
void track_mpi_all_reduce (const std::string& prefix, const ekat::Comm& comm,
                           const T* const my_vals, T* const vals, const int count, const MPI_Op op)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("mpi");
  ts.start_timer("mpi::all_reduce");
  ts.start_timer(prefix + "::mpi::all_reduce");
  comm.all_reduce(my_vals, vals, count, op);
  ts.stop_timer("mpi");
  ts.stop_timer("mpi::all_reduce");
  ts.stop_timer(prefix + "::mpi::all_reduce");
}

// Use MPI_IN_PLACE
template<typename T>
void track_mpi_all_reduce (const std::string& prefix, const ekat::Comm& comm,
                           T* const vals, const int count, const MPI_Op op)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("mpi");
  ts.start_timer("mpi::all_reduce");
  ts.start_timer(prefix + "::mpi::all_reduce");
  comm.all_reduce(vals, count, op);
  ts.stop_timer("mpi");
  ts.stop_timer("mpi::all_reduce");
  ts.stop_timer(prefix + "::mpi::all_reduce");
}

} // namespace cldera

#endif // CLDERA_MPI_TIMING_WRAPPERS_HPP
