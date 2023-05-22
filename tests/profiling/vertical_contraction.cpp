#include "profiling/stats/cldera_field_vertical_contraction.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <catch2/catch.hpp>

#include <map>

TEST_CASE ("vertical_contraction") {
  using namespace cldera;

  // REQUIRE_THROWS causes some start_timer calls to not be matched
  // by a corresponding stop_timer. When that happens, a subsequent
  // call to start_timer for that timer would throw. For this test,
  // we can just disable timings
  timing::TimingSession::instance().toggle_session(false);

  ekat::Comm comm(MPI_COMM_WORLD);

  int ymd = 20220915;
  int tod = 43000;
  TimeStamp time(ymd,tod);

  // Field dimensions
  constexpr int nlevs = 4;
  constexpr int ndims = 3;
  constexpr int ncols = 4;
  constexpr int nparts = 2;
  constexpr int part_size = ncols / nparts;

  // Create input fields
  Field s2d("scalar2d", {            ncols}, {              "ncol"}, nparts, 0,DataAccess::Copy);
  Field s3d("scalar3d", {nlevs,      ncols}, {"lev",        "ncol"}, nparts, 1,DataAccess::Copy);
  Field v3d("vector3d", {nlevs,ndims,ncols}, {"lev", "dim", "ncol"}, nparts, 2,DataAccess::Copy);
  for (int i = 0; i < nparts; ++i) {
    s2d.set_part_size(i, part_size);
    s3d.set_part_size(i, part_size);
    v3d.set_part_size(i, part_size);
  }
  s2d.commit();
  s3d.commit();
  v3d.commit();

  for (int i = 0; i < nparts; ++i) {
    std::vector<Real> data;
    // s2d
    data.resize(part_size);
    std::iota(data.begin(), data.end(), i*part_size);
    s2d.copy_part_data(i, data.data());

    // s3d
    data.resize(part_size*nlevs);
    for (int lev=0; lev<nlevs; ++lev) {
      auto beg = data.begin()+lev*part_size;
      auto end = beg + part_size;
      std::iota(beg,end,lev*ncols + i*part_size);
    }
    s3d.copy_part_data(i, data.data());

    // v3d
    data.resize(part_size*nlevs*ndims);
    for (int lev=0; lev<nlevs; ++lev) {
      for (int dim=0; dim<ndims; ++dim) {
        auto beg = data.begin()+lev*ndims*part_size + dim*part_size;
        auto end = beg + part_size;
        std::iota(beg,end,lev*ndims*ncols + dim*ncols + i*part_size);
      }
    }
    v3d.copy_part_data(i, data.data());
  }

  // Create stat
  ekat::ParameterList pl("vert_int");
  std::vector<int> lev_bounds{1,2};
  pl.set("level_bounds",lev_bounds);
  pl.set("average",false);
  auto stat = std::make_shared<FieldVerticalContraction>(comm,pl);

  // Input field does not have vertical dim
  REQUIRE_THROWS (stat->set_field(s2d));

  // Test s3d
  stat->set_field(s3d);
  auto s3d_vi = stat->compute(time);
  REQUIRE (s3d_vi.layout().rank()==1);
  auto s3d_vi_view = s3d_vi.nd_view<Real,1>();
  for (int icol=0; icol<ncols; ++icol) {
    Real tgt = 0;
    for (int ilev=lev_bounds[0];ilev<=lev_bounds[1]; ++ilev) {
      tgt += ilev*ncols+icol;
    }
    REQUIRE (s3d_vi_view(icol)==tgt);
  }

  // Test v3d
  stat->set_field(v3d);
  auto v3d_vi = stat->compute(time);
  REQUIRE (v3d_vi.layout().rank()==2);
  auto v3d_vi_view = v3d_vi.nd_view<Real,2>();
  for (int icol=0; icol<ncols; ++icol) {
    for (int idim=0; idim<ndims; ++idim) {
      Real tgt = 0;
      for (int ilev=lev_bounds[0];ilev<=lev_bounds[1]; ++ilev) {
        tgt += ilev*ncols*ndims + idim*ncols + icol;
      }
      REQUIRE (v3d_vi_view(idim,icol)==tgt);
    }
  }

}
