#ifndef CLDERA_PROFILING_TYPES_HPP
#define CLDERA_PROFILING_TYPES_HPP

#include <vector>

namespace cldera {

using Real = double;

// A field is a pointer to the array data,
// together with the info necessary to
// reshape it to an Nd array.
struct Field {
  const Real*       data;
  std::vector<int>  dims;

  int rank () const { return dims.size(); };
};

// Type of statistics to be tracked on a Field
enum class StatType {
  Max,
  Min,
  Avg,
  Sum
};

// A small struct, containing history for
// the statistics of a Field
struct History {
  StatType            stat;
  std::vector<Real>   times;
  std::vector<Real>   values;
};

} // namespace cldera

#endif // CLDERA_PROFILING_TYPES_HPP
