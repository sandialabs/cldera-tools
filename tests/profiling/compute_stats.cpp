#include "profiling/stats/cldera_compute_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <catch2/catch.hpp>

#include <map>

TEST_CASE ("stats") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  constexpr auto Max = StatType::Max;
  constexpr auto Min = StatType::Min;
  constexpr auto Avg = StatType::Avg;
  constexpr auto Sum = StatType::Sum;

  constexpr int dim0 = 2;
  constexpr int dim1 = 3;

  // Create raw data
  double raw_data[dim0*dim1];
  for (int i=0; i<dim0*dim1; ++i) {
    raw_data[i] = i;
  }

  Field f("f",{dim0,dim1});
  f.set_data(raw_data);
  std::map<StatType,History> hist;

  // Compute all supported stats
  auto stats = {Max,Min,Avg,Sum};
  for (auto stat : stats) {
    compute_stat(1.0,f,stat,hist[stat],comm);
  }

  std::map<StatType,std::vector<Real>> expected;
  expected[Max] = {5};
  expected[Min] = {0};
  expected[Sum] = {15};
  expected[Avg] = {15.0/6};

  for (auto stat : stats) {
    REQUIRE (hist[stat].times()==std::vector<Real>{1});
    REQUIRE (hist[stat].values()==expected[stat]);
  }
}
