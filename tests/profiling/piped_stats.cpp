#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/stats/cldera_register_stats.hpp"

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

  register_stats();

  ekat::Comm comm(MPI_COMM_WORLD);

  int ymd = 20220915;
  int tod = 43000;
  TimeStamp time(ymd,tod);

  // Field dimensions
  constexpr int nlevs = 4;
  constexpr int ncols = 4;
  constexpr int nparts = 2;
  constexpr int part_size = ncols / nparts;

  // Create input field
  Field s3d("scalar3d", {nlevs,      ncols}, {"lev",        "ncol"}, nparts, 1,DataAccess::Copy);
  for (int i = 0; i < nparts; ++i) {
    s3d.set_part_extent(i, part_size);
  }
  s3d.commit();

  std::vector<Real> data(part_size*nlevs);
  for (int i = 0; i < nparts; ++i) {
    for (int lev=0; lev<nlevs; ++lev) {
      auto beg = data.begin()+lev*part_size;
      auto end = beg + part_size;
      std::iota(beg,end,lev*ncols + i*part_size);
    }
    s3d.copy_part_data(i, data.data());
  }

  // Compute global_sum and compare against a piped vertical_contraction+sum_along_cols
  auto& f = StatFactory::instance();

  ekat::ParameterList pl_gsum ("gsum");
  pl_gsum.set<std::string>("type","global_sum");
  auto gsum = f.create("global_sum",comm,pl_gsum);

  ekat::ParameterList pl_pipe ("piped");
  pl_pipe.set<std::string>("type","pipe");
  pl_pipe.sublist("outer").set<std::string>("type","sum_along_columns");
  pl_pipe.sublist("inner").set<std::string>("type","vertical_contraction");
  pl_pipe.sublist("inner").set<std::vector<int>>("level_bounds",{0,nlevs-1}); // All levels
  pl_pipe.sublist("inner").set("average",false);
  auto pipe = f.create("pipe",comm,pl_pipe);

  gsum->set_field(s3d);
  gsum->create_stat_field();
  pipe->set_field(s3d);
  pipe->create_stat_field();

  auto stat_gsum = gsum->compute(time);
  auto stat_pipe = pipe->compute(time);

  REQUIRE (stat_pipe.layout()==stat_gsum.layout());
  REQUIRE (stat_gsum.data<Real>()[0] == stat_pipe.data<Real>()[0]);
}
