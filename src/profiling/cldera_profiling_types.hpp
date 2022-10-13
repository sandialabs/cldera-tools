#ifndef CLDERA_PROFILING_TYPES_HPP
#define CLDERA_PROFILING_TYPES_HPP

#include "cldera_config.h"
#include "cldera_time_stamp.hpp"
#include <ekat/ekat_assert.hpp>
#include <ekat/kokkos/ekat_kokkos_types.hpp>

#include <vector>
#include <iostream>

namespace cldera {

using Real = double;

using KokkosTypesHost = ekat::KokkosTypes<ekat::HostDevice>;

template<typename T, typename MT = Kokkos::MemoryManaged>
using view_1d_host = typename KokkosTypesHost::template view_1d<T,MT>;

// A pair used to store min/max bounds
struct Bounds {
  Real min;
  Real max;
};

// Type of statistics to be tracked on a Field
enum class StatType {
  Max,
  Min,
  Avg,
  Sum,
  Bad
};

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
