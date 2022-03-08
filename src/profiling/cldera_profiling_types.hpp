#ifndef CLDERA_PROFILING_TYPES_HPP
#define CLDERA_PROFILING_TYPES_HPP

#include <ekat/ekat_assert.hpp>

#include <vector>
#include <string>

namespace cldera {

using Real = double;

// A field is a pointer to the array data,
// together with the info necessary to
// reshape it to an Nd array.
struct Field {
  const Real*       data;
  std::vector<int>  dims;

  int rank () const { return dims.size(); };
  int size () const {
    int s = 1;
    for (auto d : dims) {
      s *= d;
    }
    return s;
  }
};

// Type of statistics to be tracked on a Field
enum class StatType {
  Max,
  Min,
  Avg,
  Sum,
  Bad
};

// A small struct, containing history for
// the statistics of a Field
// NOTE: a std;:pair would do too, but times/values are more meaningful
//       names than first/second. Besides, we can expand the struct,
//       without having to change downstream code.
struct History {
  std::vector<Real>   times;
  std::vector<Real>   values;
};

} // namespace cldera

#endif // CLDERA_PROFILING_TYPES_HPP
