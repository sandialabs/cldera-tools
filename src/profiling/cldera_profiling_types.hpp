#ifndef CLDERA_PROFILING_TYPES_HPP
#define CLDERA_PROFILING_TYPES_HPP

#include "cldera_config.h"
#include <ekat/ekat_assert.hpp>
#include <vector>

namespace cldera {

using Real = double;


// Type of statistics to be tracked on a Field
enum class StatType {
  Max,
  Min,
  Avg,
  Sum,
  Bad
};

// An std::pair would work too, but ymd/tod convey
// more meaning than first/second.
struct TimeStamp {
  int ymd;
  int tod;

};

inline bool operator== (const TimeStamp& lhs, const TimeStamp& rhs) {
  return lhs.ymd==rhs.ymd && lhs.tod==rhs.tod;
}

// A small struct, containing history for
// the statistics of a Field
// NOTE: a std;:pair would do too, but times/values are more meaningful
//       names than first/second. Besides, we can expand the struct,
//       without having to change downstream code.
class History {
public:
  int size () const { return m_times.size(); }

  void store (const int ymd, const int tod, const Real v) {
    store ({ymd,tod},v);
  }

  void store (const TimeStamp& t, const Real v) {
    m_times.push_back(t);
    m_values.push_back(v);
  }

  const std::vector<TimeStamp>& times  () const { return m_times ; }
  const std::vector<Real>&      values () const { return m_values; }

private:

  std::vector<TimeStamp>  m_times;
  std::vector<Real>       m_values;
};

} // namespace cldera

#endif // CLDERA_PROFILING_TYPES_HPP
