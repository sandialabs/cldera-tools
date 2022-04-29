#ifndef CLDERA_COMPUTE_STAT_HPP
#define CLDERA_COMPUTE_STAT_HPP

#include "profiling/stats/cldera_stats_utils.hpp"
#include "profiling/cldera_profiling_types.hpp"
#include "profiling/cldera_field.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

Real compute_max (const Field& f, const ekat::Comm& comm);
Real compute_min (const Field& f, const ekat::Comm& comm);
Real compute_sum (const Field& f, const ekat::Comm& comm);
Real compute_avg (const Field& f, const ekat::Comm& comm);

inline void compute_stat (
    const TimeStamp& t,
    const Field& f,
    const StatType s,
    History& hist,
    const ekat::Comm& comm)
{
  Real stat;
  switch (s) {
    case StatType::Max:
      stat = compute_max (f,comm); break;
    case StatType::Min:
      stat = compute_min (f,comm); break;
    case StatType::Sum:
      stat = compute_sum (f,comm); break;
    case StatType::Avg:
      stat = compute_avg (f,comm); break;
    default:
      EKAT_ERROR_MSG ("[compute_stat] Error! Stat '" + e2str(s) + "' not yet supported.\n");
  }

  // Update the history
  hist.store(t,stat);
}

} // namespace cldera

#endif // CLDERA_COMPUTE_STAT_HPP
