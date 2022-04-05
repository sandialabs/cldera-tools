#include "profiling/cldera_field.hpp"

#include <catch2/catch.hpp>

TEST_CASE ("field") {
  using namespace cldera;

  std::vector<Real> foo_data (5);
  std::vector<Real> bar_data (20);
  std::vector<std::vector<Real>> baz_data (3,std::vector<Real>(20));

  std::iota(foo_data.begin(),foo_data.end(),0);
  std::iota(bar_data.begin(),bar_data.end(),0);
  std::iota(baz_data[0].begin(),baz_data[0].end(),0);
  std::iota(baz_data[1].begin(),baz_data[1].end(),20);
  std::iota(baz_data[2].begin(),baz_data[2].end(),40);

  Field foo("foo",{5},foo_data.data());
  Field bar("bar",{5,4});
  Field baz("baz",{5,4,3},3,2); // 3 partitions along last dim
  REQUIRE_THROWS (Field("",{5},6,0)); // Too many parts
  REQUIRE_THROWS (Field("",{5},1,1)); // part dim OOB
  REQUIRE (foo.committed());

  // Set data in fields
  bar.set_data(bar_data.data());
  baz.set_part_data(0,0,baz_data[0].data());
  baz.set_part_data(2,2,baz_data[2].data());
  baz.set_part_data(1,1,baz_data[1].data());
  REQUIRE_THROWS (foo.set_data(foo_data.data())); // Data already set
  REQUIRE_THROWS (baz.set_part_data(3,0,foo_data.data())); // part idx OOB
  REQUIRE_THROWS (baz.set_part_data(1,10,foo_data.data())); // part beg OOB

  // Get data
  REQUIRE_THROWS (baz.get_part_data(0)); // Field not yet committed
  bar.commit();
  baz.commit();
  REQUIRE_THROWS (baz.get_data()); // Can't call get_data on partitioned field

  // Check dimensions
  REQUIRE(foo.layout().size()==5);
  REQUIRE(bar.layout().size()==20);
  REQUIRE(baz.layout().size()==60);

  REQUIRE(baz.part_layout(0).size()==20);
  REQUIRE(baz.part_layout(1).size()==20);
  REQUIRE(baz.part_layout(2).size()==20);

  // Check data
  for (int i=0; i<foo.layout().size(); ++i) {
    REQUIRE (foo.get_data()[i]==i);
  }
  for (int i=0; i<bar.layout().size(); ++i) {
    REQUIRE (bar.get_data()[i]==i);
  }
  for (int p=0,offset=0; p<baz.nparts(); ++p) {
    for (int i=0; i<baz.part_layout(p).size(); ++i) {
      REQUIRE (baz.get_part_data(p)[i]==(i+offset));
    }
    offset += baz.part_layout(p).size();
  }
}
