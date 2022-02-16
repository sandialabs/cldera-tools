#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_archive.hpp"
#include "Kokkos_Core.hpp"

template<typename DT>
using View = Kokkos::View<DT,Kokkos::HostSpace>;

TEST_CASE ("archive") {
  using namespace cldera;

  ProfilingArchive archive;

  View<Real>    foo ("foo");
  View<Real*>   bar ("bar",5);
  View<Real**>  baz ("baz",5,4);
  View<Real***> qux ("qux",5,4,3);

  // Invalid pointer
  REQUIRE_THROWS (archive.add_field("foo",nullptr,{1}));

  archive.add_field("foo",foo.data(),{});
  archive.add_field("bar",bar.data(),{5});
  archive.add_field("baz",baz.data(),{5,4});
  archive.add_field("qux",qux.data(),{5,4,3});

  // Field already registered
  REQUIRE_THROWS (archive.add_field("qux",qux.data(),{5,4,3}));

  // Check pointers are set correctly
  foo() = 10;
  REQUIRE (archive.get_field("foo").data[0]==10);

  // Store a bunch of data
  constexpr auto Avg = StatType::Avg;
  constexpr auto Max = StatType::Max;
  constexpr auto Min = StatType::Min;
  constexpr auto Sum = StatType::Sum;
  for (Real t=0; t<5.0; t+=1.0) {
    for (auto stat : {Avg,Max,Min,Sum}) {
      archive.store_stat("foo",stat,t,2*t);
    }
  }
  // Time must increase
  REQUIRE_THROWS(archive.store_stat("foo",StatType::Max,0,0));
  // Field not registered
  REQUIRE_THROWS(archive.store_stat("xyz",StatType::Max,0,0));

  // Check history data
  std::vector<Real> t = {0,1,2,3,4};
  std::vector<Real> v = {0,2,4,6,8};
  for (auto stat : {Avg,Max,Min,Sum}) {
    const auto& history = archive.get_stat_history("foo",stat);
    REQUIRE (history.times.size()==5);
    REQUIRE (history.values.size()==5);

    REQUIRE (history.times==t);
    REQUIRE (history.values==v);
  }
}
