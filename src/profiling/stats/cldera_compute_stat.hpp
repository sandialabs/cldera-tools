#ifndef CLDERA_COMPUTE_STAT_HPP
#define CLDERA_COMPUTE_STAT_HPP

#include "profiling/stats/cldera_stats_utils.hpp"
#include "profiling/cldera_profiling_types.hpp"
#include "profiling/cldera_field.hpp"

namespace cldera {

void compute_max (const Field& f, History& hist);
void compute_min (const Field& f, History& hist);
void compute_sum (const Field& f, History& hist);
void compute_avg (const Field& f, History& hist);

inline void compute_stat (const Real time, const Field& f, const StatType s, History& hist)
{
  switch (s) {
    case StatType::Max:
      compute_max (f,hist); break;
    case StatType::Min:
      compute_min (f,hist); break;
    case StatType::Sum:
      compute_sum (f,hist); break;
    case StatType::Avg:
      compute_avg (f,hist); break;
    default:
      EKAT_ERROR_MSG ("[compute_stat] Error! Stat '" + e2str(s) + "' not yet supported.\n");
  }
  // Update the times
  hist.times.push_back(time);
}

} // namespace cldera

#endif // CLDERA_COMPUTE_STAT_HPP
