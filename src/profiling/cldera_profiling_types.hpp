#ifndef CLDERA_PROFILING_TYPES_HPP
#define CLDERA_PROFILING_TYPES_HPP

#include "cldera_config.h"
#include "cldera_time_stamp.hpp"
#include <ekat/ekat_assert.hpp>
#include <ekat/kokkos/ekat_kokkos_types.hpp>
#include <ekat/std_meta/ekat_std_type_traits.hpp>

#include <vector>
#include <iostream>

namespace ekat {
// EKAT's DataND recursion has N=1 as base case, causing infinite recursion if N=0
// TODO: fix in EKAT
template<typename T>
struct DataND<T,0> {
  using type = T;
};
} // namespace ekat

namespace cldera {

using Real = double;

using KokkosTypesHost = ekat::KokkosTypes<ekat::HostDevice>;

template<typename T, typename MT = Kokkos::MemoryManaged>
using view_1d_host = typename KokkosTypesHost::template view_1d<T,MT>;

template<typename T, int N, typename MT = Kokkos::MemoryManaged>
using view_Nd_host = typename KokkosTypesHost::template view_ND<T,N,MT>;

// A pair used to store min/max bounds
template<typename T>
struct Bounds {
  Bounds () = default;
  Bounds (const T min_, const T max_)
   : min(min_), max(max_) {}

  Bounds (const std::vector<T>& bounds) {
    EKAT_REQUIRE_MSG (bounds.size()==2,
        "Error! Invalid size for bounds vector (" << bounds.size() << ").\n");
    min = bounds[0];
    max = bounds[1];
  }

  std::vector<T> to_vector () const {
    return std::vector<T>{min,max};
  }

  bool contains (const T val) {
    return val>=min && val<=max;
  }

  bool contains (const T val, bool closed_left, bool closed_right) {
    return (closed_left  ? val>=min : val>min) and
           (closed_right ? val<=max : val<max);
  }

  T min;
  T max;
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
