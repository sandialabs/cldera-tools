#include "profiling/cldera_profiling_types.hpp"

#include <limits>

namespace cldera {

void compute_min (const Field& f, History& hist)
{
  EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
      "Error! The type cldera::Real is not capable of representing infinity.\n");

  Real min = std::numeric_limits<Real>::infinity();
  for (int i=0; i<f.size(); ++i) {
    min = std::min(min,f.data[i]);
  }

  hist.values.push_back(min);
}

} // namespace cldera
