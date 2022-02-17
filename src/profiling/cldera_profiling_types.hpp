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
};

// Type of statistics to be tracked on a Field
enum class StatType {
  Max,
  Min,
  Avg,
  Sum,
  Bad
};

inline std::string e2str (const StatType s)
{
  switch (s) {
    case StatType::Max: return "Max";
    case StatType::Min: return "Min";
    case StatType::Avg: return "Avg";
    case StatType::Sum: return "Sum";
    case StatType::Bad: return "INVALID";
    default:
      // Note: this should not happen, but if we add enum values, without
      //       adding a 'case' to this function, AND without paying attention
      //       to compiler warnings, we'd get long segfault messages.
      //       This way, we at least get an intelligible message.
      EKAT_ERROR_MSG ("Unhandled StatType value.\n");
  }
}

inline StatType str2stat (const std::string& str)
{
  for (auto s : {StatType::Max, StatType::Min, StatType::Avg, StatType::Sum}) {
    if (str==e2str(s)) {
      return s;
    }
  }

  EKAT_ERROR_MSG ("Unrecognized string for conversion to StatType enum.\n");

  return StatType::Bad;
}


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
