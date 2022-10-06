#ifndef CLDERA_PROFILING_TYPES_HPP
#define CLDERA_PROFILING_TYPES_HPP

#include "cldera_config.h"
#include <ekat/ekat_assert.hpp>
#include <ekat/kokkos/ekat_kokkos_types.hpp>

#include <vector>
#include <iostream>

namespace cldera {

using Real = double;

using KokkosTypesHost = ekat::KokkosTypes<ekat::HostDevice>;

template<typename T, typename MT = Kokkos::MemoryManaged>
using view_1d_host = typename KokkosTypesHost::template view_1d<T,MT>;

} // namespace cldera

#endif // CLDERA_PROFILING_TYPES_HPP
