#include "profiling/cldera_field.hpp"

#include <catch2/catch.hpp>

TEST_CASE ("field") {
  using namespace cldera;

  std::vector<Real> foo_data (5);
  std::vector<int> bar_data (20);
  std::vector<std::vector<Real>> baz_data (3,std::vector<Real>(20));

  std::iota(foo_data.begin(),foo_data.end(),0.0);
  std::iota(bar_data.begin(),bar_data.end(),0);
  std::iota(baz_data[0].begin(),baz_data[0].end(),0.0);
  std::iota(baz_data[1].begin(),baz_data[1].end(),20.0);
  std::iota(baz_data[2].begin(),baz_data[2].end(),40.0);

  Field foo("foo",{5},{"col"},foo_data.data());
  Field bar("bar",{5,4},{"col","lev"},DataAccess::View,DataType::IntType);
  Field baz("baz",{5,4,3},{"col","cmp","lev"},3,2); // 3 partitions along last dim
  Field foobar("foobar",{5},{"col"},DataAccess::Copy);
  REQUIRE_THROWS (Field("",{5},{"col"},6,0)); // Too many parts
  REQUIRE_THROWS (Field("",{5},{"col"},1,1)); // part dim OOB
  REQUIRE (foo.committed());

  // Test FieldLayout::to_string()
  REQUIRE (foo.layout().to_string()=="<col> (5)");
  REQUIRE (bar.layout().to_string()=="<col,lev> (5,4)");
  REQUIRE (baz.layout().to_string()=="<col,cmp,lev> (5,4,3)");

  foobar.set_part_size(0,5); // OK, same value
  foobar.set_part_size(0,5); // OK, same value
  REQUIRE_THROWS(foobar.set_part_size(0,6)); // Can't change part sizes

  // Set data in fields
  REQUIRE_THROWS (foo.set_data(foo_data.data())); // Data already set
  REQUIRE_THROWS (baz.set_part_size(3,1)); // part idx OOB
  REQUIRE_THROWS (baz.set_part_size(1,10)); // part size OOB
  REQUIRE_THROWS (foobar.set_part_data<Real>(0,foo_data.data())); // Not a View field
  REQUIRE_THROWS (bar.set_data<int>(nullptr)); // Invalid pointer
  REQUIRE_THROWS (bar.set_data(foo_data.data())); // Invalid data type
  bar.set_data(bar_data.data());
  baz.set_part_size(0,1);
  baz.set_part_size(2,1);
  baz.set_part_size(1,1);
  baz.set_part_data(0,baz_data[0].data());
  baz.set_part_data(2,baz_data[2].data());
  baz.set_part_data(1,baz_data[1].data());

  // Commit
  REQUIRE_THROWS (baz.part_data<Real>(0)); // Field not yet committed
  bar.commit();
  baz.commit();
  foobar.commit();
  REQUIRE_THROWS (bar.data<Real>()); // Wrong data type
  REQUIRE_THROWS (bar.set_data(foo_data.data())); // Can't reset data pointer
  REQUIRE_THROWS (baz.data<Real>()); // Can't call data on partitioned field

  // Copy data
  REQUIRE_THROWS (foo.copy_data(foo_data.data())); // Can't copy in a View field.
  foobar.copy_data(foo_data.data());

  // Check dimensions
  REQUIRE(foo.layout().size()==5);
  REQUIRE(bar.layout().size()==20);
  REQUIRE(baz.layout().size()==60);

  REQUIRE(baz.part_layout(0).size()==20);
  REQUIRE(baz.part_layout(1).size()==20);
  REQUIRE(baz.part_layout(2).size()==20);

  // Check data
  for (int i=0; i<foo.layout().size(); ++i) {
    REQUIRE (foo.data<Real>()[i]==i);
    REQUIRE (foo.data<Real>()[i]==foobar.data<Real>()[i]);
  }
  for (int i=0; i<bar.layout().size(); ++i) {
    REQUIRE (bar.data<int>()[i]==i);
  }
  for (int p=0,offset=0; p<baz.nparts(); ++p) {
    for (int i=0; i<baz.part_layout(p).size(); ++i) {
      REQUIRE (baz.part_data<Real>(p)[i]==(i+offset));
    }
    offset += baz.part_layout(p).size();
  }

  // Check degenerate case of a scalar
  Field scalar("s",{},{},DataAccess::Copy);

  REQUIRE(scalar.layout().rank()==0);
  REQUIRE(scalar.layout().size()==1);
  REQUIRE(scalar.part_layout(0).rank()==0);
  REQUIRE(scalar.part_layout(0).size()==1);

  // Check clone method
  auto baz2 = baz.clone();
  REQUIRE (baz2.name()==baz.name());
  REQUIRE (baz2.layout()==baz.layout());
  REQUIRE (baz2.nparts()==baz.nparts());
  REQUIRE (baz2.part_dim()==baz.part_dim());
  REQUIRE (baz2.data_access()==DataAccess::Copy);
  REQUIRE (baz2.data_type()==baz.data_type());
  for (int p=0; p<baz.nparts(); ++p) {
    REQUIRE (baz2.part_layout(p)==baz.part_layout(p));
    REQUIRE (baz2.part_layout(p)==baz.part_layout(p));

    auto baz_pv  = baz.part_view<Real>(p);
    auto baz2_pv = baz2.part_view<Real>(p);
    REQUIRE (baz_pv.size()==baz2_pv.size());
    REQUIRE (baz_pv.data()!=baz2_pv.data()); // Check it's not a shallow clone!
    for (size_t i=0; i<baz_pv.size(); ++i) {
      REQUIRE (baz_pv[i]==baz2_pv[i]);
    }
  }

  // Check deep copy method
  baz2.deep_copy(-1.0);
  baz2.deep_copy(baz);
  for (int p=0; p<baz.nparts(); ++p) {
    auto baz_pv  = baz.part_view<Real>(p);
    auto baz2_pv = baz2.part_view<Real>(p);
    for (size_t i=0; i<baz_pv.size(); ++i) {
      REQUIRE (baz_pv[i]==baz2_pv[i]);
    }
  }
}
