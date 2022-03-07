#ifndef CLDERA_STATS_UTILS_HPP
#define CLDERA_STATS_UTILS_HPP

#include "profiling/cldera_profiling_types.hpp"

namespace cldera {

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

} // namespace cldera

#endif // CLDERA_STATS_UTILS_HPP
