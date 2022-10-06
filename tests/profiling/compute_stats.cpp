#include "profiling/stats/cldera_field_stat_factory.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <catch2/catch.hpp>

#include <map>

TEST_CASE ("stats") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  constexpr int dim0 = 2;
  constexpr int dim1 = 3;

  // Create raw data
  double raw_data[dim0*dim1];
  for (int i=0; i<dim0*dim1; ++i) {
    raw_data[i] = i;
  }

  Field f("f",{dim0,dim1},{"dim0","dim1"});
  f.set_data(raw_data);

  // Compute all supported stats
  auto stats_names = {"global_max", "global_min", "global_sum", "global_avg"};
  std::map<std::string,Field> stat_fields;
  std::map<std::string,Field> expected;
  for (auto sname : stats_names) {
    auto stat = create_stat(sname,comm);
    stat_fields.emplace(sname,stat->compute(f));

    expected.emplace(sname,Field(sname,stat_fields.at(sname).layout(),DataAccess::Copy));
  }

  expected["global_max"].data_nonconst()[0] = 5;
  expected["global_min"].data_nonconst()[0] = 0;
  expected["global_sum"].data_nonconst()[0] = 15;
  expected["global_avg"].data_nonconst()[0] = 15.0/6;

  for (auto sname : stats_names) {
    REQUIRE (expected.at(sname).data()[0]==stat_fields.at(sname).data()[0]);
  }
}
