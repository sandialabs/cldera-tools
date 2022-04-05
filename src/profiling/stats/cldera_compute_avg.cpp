#include "profiling/stats/cldera_compute_stat.hpp"
#include "profiling/cldera_field.hpp"

namespace cldera {

void compute_avg (const Field& f, History& hist)
{
  compute_sum (f,hist);
  hist.values.back() /= f.layout().size();
}

} // namespace cldera
