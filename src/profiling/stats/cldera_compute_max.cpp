#include "profiling/cldera_profiling_types.hpp"

#include <limits>

namespace cldera {

void compute_max (const Field& f, History& hist)
{
  EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
      "Error! The type cldera::Real is not capable of representing infinity.\n");

  Real max = -std::numeric_limits<Real>::infinity();
  for (int i=0; i<f.size(); ++i) {
    max = std::max(max,f.data[i]);
  }

  hist.values.push_back(max);
}

} // namespace cldera
