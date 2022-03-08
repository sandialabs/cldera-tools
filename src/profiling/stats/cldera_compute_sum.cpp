#include "profiling/cldera_profiling_types.hpp"
#include "profiling/stats/cldera_compute_stat.hpp"

namespace cldera {

void compute_sum (const Field& f, History& hist)
{
  // Note: use Kahan summation to increase accuracy
  Real sum = 0;
  Real c = 0;
  Real temp,y;
  for (int i=0; i<f.size(); ++i) {
    y = f.data[i] - c;
    temp = sum + y;
    c = (temp - sum) - y;
    sum = temp;
  }
  hist.values.push_back(sum);
}

} // namespace cldera
