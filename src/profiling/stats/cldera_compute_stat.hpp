#ifndef CLDERA_COMPUTE_STAT_HPP
#define CLDERA_COMPUTE_STAT_HPP

#include "profiling/stats/cldera_stats_utils.hpp"
#include "profiling/cldera_profiling_types.hpp"
#include "profiling/cldera_field.hpp"

namespace cldera {

Real compute_max (const Field& f);
Real compute_min (const Field& f);
Real compute_sum (const Field& f);
Real compute_avg (const Field& f);

inline void compute_stat (const Real time, const Field& f, const StatType s, History& hist)
{
  Real stat;
  switch (s) {
    case StatType::Max:
      stat = compute_max (f); break;
    case StatType::Min:
      stat = compute_min (f); break;
    case StatType::Sum:
      stat = compute_sum (f); break;
    case StatType::Avg:
      stat = compute_avg (f); break;
    default:
      EKAT_ERROR_MSG ("[compute_stat] Error! Stat '" + e2str(s) + "' not yet supported.\n");
  }

  // Update the history
  hist.store(time,stat);
}

} // namespace cldera

#endif // CLDERA_COMPUTE_STAT_HPP
