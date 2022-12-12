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
    expected.at(sname).commit();
  }

  expected["global_max"].data_nonconst<Real>()[0] = 5;
  expected["global_min"].data_nonconst<Real>()[0] = 0;
  expected["global_sum"].data_nonconst<Real>()[0] = 15;
  expected["global_avg"].data_nonconst<Real>()[0] = 15.0/6;

  for (auto sname : stats_names) {
    REQUIRE (expected.at(sname).data<Real>()[0]==stat_fields.at(sname).data<Real>()[0]);
  }
}

TEST_CASE ("stats along columns") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  constexpr int dim0 = 2;
  constexpr int dim1 = 3;
  constexpr int dim2 = 4;

  // Create raw data
  double raw_data[dim0*dim1*dim2];
  for (int i=0; i<dim0*dim1*dim2; ++i) {
    raw_data[i] = i;
  }

  Field f("f", {dim0,dim1,dim2}, {"lev", "dim", "ncol"});
  f.set_data(raw_data);

  auto stats_names = {"max_along_columns", "min_along_columns",
      "sum_along_columns", "avg_along_columns"};
  std::map<std::string, Field> stat_fields;
  std::map<std::string, Field> expected;
  for (auto sname : stats_names) {
    auto stat = create_stat(sname, comm);
    stat_fields.emplace(sname, stat->compute(f));
    expected.emplace(sname, Field(sname, stat_fields.at(sname).layout(), DataAccess::Copy));
    expected.at(sname).commit();
  }

  Real max_values[] = {3,7,11,15,19,23};
  Real min_values[] = {0,4,8,12,16,20};
  Real sum_values[] = {6,22,38,54,70,86};
  Real avg_values[] = {6.0/4,22.0/4,38.0/4,54.0/4,70.0/4,86.0/4};
  for (int i = 0; i < dim0*dim1; ++i) {
    expected["max_along_columns"].data_nonconst<Real>()[i] = max_values[i];
    expected["min_along_columns"].data_nonconst<Real>()[i] = min_values[i];
    expected["sum_along_columns"].data_nonconst<Real>()[i] = sum_values[i];
    expected["avg_along_columns"].data_nonconst<Real>()[i] = avg_values[i];
  }

  for (auto sname : stats_names)
    for (int i = 0; i < dim0*dim1; ++i)
      REQUIRE (expected.at(sname).data<Real>()[i]==stat_fields.at(sname).data<Real>()[i]);
}

TEST_CASE ("stats along columns with parts") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  // Allocate field with multiple parts
  constexpr int dim0 = 2;
  constexpr int dim1 = 3;
  constexpr int dim2 = 4;
  constexpr int nparts = 4;
  constexpr int part_dim = 2;
  Field foo("foo", {dim0,dim1,dim2}, {"lev", "dim", "ncol"}, nparts, part_dim);
  std::vector<std::vector<Real>> foo_data (dim2,std::vector<Real>(dim0*dim1));
  for (int i = 0; i < dim2; ++i) {
    std::iota(foo_data[i].begin(), foo_data[i].end(), i*dim0*dim1);
    const int part_size = 1;
    foo.set_part_size(i, part_size);
    foo.set_part_data(i, foo_data[i].data());
  }
  foo.commit();

  // Compute all supported stats along columns
  auto stats_names = {"max_along_columns", "min_along_columns",
      "sum_along_columns", "avg_along_columns"};
  std::map<std::string, Field> stat_fields;
  std::map<std::string, Field> expected;
  for (auto sname : stats_names) {
    auto stat = create_stat(sname, comm);
    stat_fields.emplace(sname, stat->compute(foo));
    expected.emplace(sname, Field(sname, stat_fields.at(sname).layout(), DataAccess::Copy));
    expected.at(sname).commit();
  }

  Real max_values[] = {18,19,20,21,22,23};
  Real min_values[] = {0,1,2,3,4,5};
  Real sum_values[] = {36,40,44,48,52,56};
  Real avg_values[] = {36.0/4,40.0/4,44.0/4,48.0/4,52.0/4,56.0/4};
  for (int i = 0; i < dim0*dim1; ++i) {
    expected["max_along_columns"].data_nonconst<Real>()[i] = max_values[i];
    expected["min_along_columns"].data_nonconst<Real>()[i] = min_values[i];
    expected["sum_along_columns"].data_nonconst<Real>()[i] = sum_values[i];
    expected["avg_along_columns"].data_nonconst<Real>()[i] = avg_values[i];
  }

  for (auto sname : stats_names)
    for (int i = 0; i < dim0*dim1; ++i)
      REQUIRE (expected.at(sname).data<Real>()[i]==stat_fields.at(sname).data<Real>()[i]);
}
