#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_types.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

#include "ekat/kokkos/ekat_kokkos_meta.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_type_traits.hpp"
#include "ekat/std_meta/ekat_std_type_traits.hpp"

#include <random>

namespace {

TEST_CASE("slice") {
  using namespace cldera;
  using kt = KokkosTypesHost;
  kt::view_ND<double,3,Kokkos::MemoryManaged> v("",2,3,4);

  std::random_device rdev;
  const int catchRngSeed = Catch::rngSeed();
  int seed = catchRngSeed==0 ? rdev() : catchRngSeed;
  std::mt19937_64 engine(seed);
  std::uniform_real_distribution<double> pdf(0,1);

  ekat::genRandArray(v,engine,pdf);

  auto sv1 = slice(v,0,1);
  auto sv2 = slice(v,1,1);
  auto sv3 = slice(v,2,1);

  // Check subview along each dim
  for (int i=0; i<2; ++i) {
    for (int j=0; j<3; ++j) {
      for (int k=0; k<4; ++k) {
        REQUIRE (sv1(j,k)==v(1,j,k));
        REQUIRE (sv2(i,k)==v(i,1,k));
        REQUIRE (sv3(i,j)==v(i,j,1));
      }
    }
  }

  // Check double subview commutative property
  auto sv1_3 = slice(sv1,1,1);
  auto sv3_1 = slice(sv3,0,1);
  for (int j=0; j<3; ++j) {
    REQUIRE (sv1_3(j)==sv3_1(j));
  }
}

} // empty namespace
