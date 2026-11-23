#include "profiling/cldera_field.hpp"

#include <catch2/catch.hpp>

TEST_CASE ("field") {
  using namespace cldera;

  SECTION ("basic") {
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

    foobar.set_part_extent(0,5); // OK, same value
    foobar.set_part_extent(0,5); // OK, same value
    REQUIRE_THROWS(foobar.set_part_extent(0,6)); // Can't change part sizes

    // Set data in fields
    REQUIRE_THROWS (foo.set_data(foo_data.data())); // Data already set
    REQUIRE_THROWS (baz.set_part_extent(3,1)); // part idx OOB
    REQUIRE_THROWS (baz.set_part_extent(1,10)); // part size OOB
    REQUIRE_THROWS (foobar.set_part_data<Real>(0,foo_data.data())); // Not a View field
    REQUIRE_THROWS (bar.set_data<int>(nullptr)); // Invalid pointer
    REQUIRE_THROWS (bar.set_data(foo_data.data())); // Invalid data type
    bar.set_data(bar_data.data());
    baz.set_part_extent(0,1);
    baz.set_part_extent(2,1);
    baz.set_part_extent(1,1);
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

  // Check uniform alloc size 
  SECTION ("part_alloc_size") {
    constexpr int nparts = 2;
    constexpr int ncols = 3;
    constexpr int nlevs = 2;
    constexpr int part_dim = 1; 

    // On each part, we'll fill only the first $part_size cols
    Real f_data [nparts][nlevs][ncols];

    FieldLayout fl({nlevs,ncols},{"nlevs","ncols"});
    Field f("test_alloc_size",fl,nparts,part_dim,DataAccess::View,DataType::RealType,ncols);
    f.set_part_extent(0,2);
    f.set_part_extent(1,1);

    f.set_part_data(0,&f_data[0][0][0]);
    f.set_part_data(1,&f_data[1][0][0]);
    f.commit();

    auto part = [&] (int icol) {
      return icol<2 ? 0 : 1;
    };

    int part_offset[2] = {0, 2};
    std::fill_n(&f_data[0][0][0],nparts*nlevs*ncols,-1);
    for (int icol=0,val=0; icol<ncols; ++icol) {
      int p = part(icol);
      for (int ilev=0; ilev<nlevs; ++ilev,++val) {
        f_data[p][ilev][icol-part_offset[p]] = val;
      }
    }
    
    for (int ipart=0; ipart<nparts; ++ipart) {
      auto fp = f.part_nd_view<Real,2>(ipart);
      REQUIRE (fp.extent(0)==nlevs);
      REQUIRE (fp.extent(1)==ncols);
      auto pl = f.part_layout(ipart);
      auto offset = f.part_offset(ipart);
      for (int icol=0; icol<pl.dims()[1]; ++icol) {
        for (int ilev=0; ilev<nlevs; ++ilev) {
          REQUIRE (fp(ilev,icol)==f_data[ipart][ilev][icol]);
        }
      }
    }
  }

  SECTION ("update") {
    constexpr int ncols = 3;
    constexpr int nlevs = 4;
    constexpr int nparts = 2;

    FieldLayout fl ({nlevs,ncols},{"nlevs","ncols"});
    FieldLayout fl_bad ({nlevs+1,ncols},{"nlevs","ncols"});
    Field yr("yr",fl,nparts,1,  DataAccess::Copy,DataType::RealType);
    Field yi("yi",fl,nparts,1,  DataAccess::Copy,DataType::IntType);
    Field xr("xr",fl,nparts,1,  DataAccess::Copy,DataType::RealType);
    Field b1("bad1",fl_bad,nparts,1,DataAccess::Copy,DataType::RealType);
    Field b2("bad1",fl,DataAccess::Copy,DataType::RealType);

    yr.set_part_extent(0,1);
    yi.set_part_extent(0,1);
    xr.set_part_extent(0,1);
    b1.set_part_extent(0,1);

    yr.set_part_extent(1,2);
    yi.set_part_extent(1,2);
    xr.set_part_extent(1,2);
    b1.set_part_extent(1,2);

    yr.commit();
    yi.commit();
    xr.commit();
    b1.commit();
    b2.commit();

    REQUIRE_THROWS(yr.update(b1,1,1)); // layout mismatch
    REQUIRE_THROWS(yr.update(b2,1,1)); // partitioning mismatch
    REQUIRE_THROWS(yi.update(xr,1,1)); // incompatible X-Y data types
    REQUIRE_THROWS(yi.update(yi,1.0,1.0)); // incompatible X-coeffs data types

    for (int p=0, k=0; p<nparts; ++p) {
      auto xpv = xr.part_view_nonconst<Real>(p);
      auto s = xr.part_layout(p).size();
      std::iota(xpv.data(),xpv.data()+s,Real(k));
      k += s;
    }
    yr.deep_copy(0);

    // Check that f1==scale*f2
    auto check = [](const Field& f1, const Field& f2, const int scale) {
      for (int p=0; p<nparts; ++p) {
        auto pv1 = f1.part_view<Real>(p);
        auto pv2 = f2.part_view<Real>(p);
        auto s = f1.part_layout(p).size();
        for (int i=0; i<s; ++i) {
          if ( pv1(i)!= pv2(i)*scale ) {
            return false;
          }
        }
      }
      return true;
    };

    // y = 0 + x;
    yr.update(xr,1,1);
    REQUIRE (check(yr,xr,1));

    // y = 0 + 2*x
    yr.update(xr,2,0);
    REQUIRE (check(yr,xr,2));

    // y = y + 1*x
    yr.update(xr,1,1);
    REQUIRE (check(yr,xr,3));
  }

  // BMW v2: ownership (Copy/View) vs writability (read-only/writable)
  SECTION ("writable") {
    constexpr int ncols = 3;
    constexpr int nlevs = 2;
    constexpr int nparts = 2;
    FieldLayout fl ({nlevs,ncols},{"lev","ncol"});

    // 1. View fields are read-only by default, even with non-const pointers
    std::vector<Real> d0 (nlevs*ncols,1.0), d1(nlevs*ncols,2.0);
    Field ro ("ro",fl,nparts,1,DataAccess::View,DataType::RealType,ncols);
    REQUIRE (ro.is_read_only());
    ro.set_part_extent(0,2);
    ro.set_part_extent(1,1);
    ro.set_part_data(0,d0.data());
    ro.set_part_data(1,d1.data());
    ro.commit();
    REQUIRE_THROWS (ro.part_view_nonconst<Real>(0));
    REQUIRE_THROWS (ro.deep_copy(0.0));
    REQUIRE (ro.part_view<Real>(1)(0)==2.0); // const access OK

    // 2. Writable View field: writes go straight into host-app memory
    Field rw ("rw",fl,nparts,1,DataAccess::View,DataType::RealType,ncols);
    rw.set_read_only(false);
    REQUIRE (not rw.is_read_only());
    rw.set_part_extent(0,2);
    rw.set_part_extent(1,1);
    rw.set_part_data(0,d0.data());
    rw.set_part_data(1,d1.data());
    rw.commit();
    REQUIRE_THROWS (rw.set_read_only(true)); // can't change after commit
    rw.part_nd_view_nonconst<Real,2>(1)(1,0) = 42.0;
    REQUIRE (d1[1*ncols+0]==42.0);
    rw.deep_copy(0.0);
    for (auto v : d0) REQUIRE (v==0.0);
    for (auto v : d1) REQUIRE (v==0.0);

    // 3. Writable View field with a const pointer: commit must fail
    Field bad ("bad",fl,DataAccess::View,DataType::RealType);
    bad.set_read_only(false);
    const Real* cptr = d0.data();
    bad.set_part_data(0,cptr);
    REQUIRE_THROWS (bad.commit());

    // 4. set_read_only must be called before setting data
    Field late ("late",fl,DataAccess::View,DataType::RealType);
    late.set_part_data(0,d0.data());
    REQUIRE_THROWS (late.set_read_only(false));

    // 5. read_only() shallow copy of a writable field refuses writes, but the
    //    original is still writable. Clones are always writable.
    auto rw_ro = rw.read_only();
    REQUIRE (rw_ro.is_read_only());
    REQUIRE_THROWS (rw_ro.part_view_nonconst<Real>(0));
    REQUIRE_NOTHROW (rw.part_view_nonconst<Real>(0));
    auto ro_clone = ro.clone();
    REQUIRE (not ro_clone.is_read_only());
    REQUIRE_NOTHROW (ro_clone.deep_copy(3.0));
    REQUIRE (d0[0]==0.0); // clone is a deep copy

    // 6. Copy fields are writable by default
    Field cp ("cp",fl,DataAccess::Copy,DataType::RealType);
    REQUIRE (not cp.is_read_only());
    cp.commit();
    REQUIRE_NOTHROW (cp.deep_copy(1.0));
  }
}
