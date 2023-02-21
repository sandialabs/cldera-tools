#include "profiling/stats/cldera_field_stat_factory.hpp"
#include "profiling/stats/cldera_field_bounding_box.hpp"
#include "profiling/stats/cldera_field_zonal_mean.hpp"

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

TEST_CASE ("stat_int_field") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  constexpr int dim0 = 10;

  // Create raw data
  int raw_data[dim0];
  for (int i=0; i<dim0; ++i) {
    raw_data[i] = i+1;
  }

  Field f("f",{dim0},{"dim0"},DataAccess::View,IntType);
  f.set_data(raw_data);

  auto stat = create_stat("global_sum",comm);
  const auto computed = stat->compute(f).data<int>()[0];
  const auto expected = (dim0+1)*dim0/2;
  REQUIRE (computed==expected);
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
  constexpr int nparts = 2;
  constexpr int part_size = dim2 / nparts;
  constexpr int part_dim = 2;
  Field foo("foo", {dim0,dim1,dim2}, {"lev", "dim", "ncol"}, nparts, part_dim);
  std::vector<std::vector<Real>> foo_data (nparts,std::vector<Real>(dim0*dim1*part_size));
  for (int i = 0; i < nparts; ++i) {
    std::iota(foo_data[i].begin(), foo_data[i].end(), i*dim0*dim1*part_size);
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

  Real max_values[] = {13,15,17,19,21,23};
  Real min_values[] = {0,2,4,6,8,10};
  Real sum_values[] = {26,34,42,50,58,66};
  Real avg_values[] = {26.0/4,34.0/4,42.0/4,50.0/4,58.0/4,66.0/4};
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

TEST_CASE ("stats - bounds") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  // Allocate field
  constexpr int dim0 = 2;
  constexpr int dim1 = 3;
  constexpr int dim2 = 4;
  constexpr int nparts = 2;
  constexpr int part_size = dim2 / nparts;
  constexpr int part_dim = 2;
  Field foo("foo", {dim0,dim1,dim2}, {"lev", "dim", "ncol"}, nparts, part_dim);
  std::vector<std::vector<Real>> foo_data (nparts,std::vector<Real>(dim0*dim1*part_size));
  for (int i = 0; i < nparts; ++i) {
    std::iota(foo_data[i].begin(), foo_data[i].end(), i*dim0*dim1*part_size);
    foo.set_part_size(i, part_size);
    foo.set_part_data(i, foo_data[i].data());
  }
  foo.commit();

  // Test bounded field
  // bounds set to 12.1, 20.1
  const auto bounded_sname = "bounded";
  const auto bounded_stat = create_stat(bounded_sname, comm);
  const auto bounded_field = bounded_stat->compute(foo);
  const Real bounded_expected[] = {
      0.0, 0.0, 0.0, 13.0,
      0.0, 0.0, 14.0, 15.0,
      0.0, 0.0, 16.0, 17.0,
      0.0, 0.0, 18.0, 19.0,
      0.0, 0.0, 20.0, 0.0,
      0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < dim0*dim1*dim2; ++i)
    REQUIRE (bounded_expected[i]==bounded_field.data<Real>()[i]);

  // Allocate lat/lon
  std::shared_ptr<Field> lat(new Field("lat", {dim2}, {"ncol"}, nparts, 0));
  std::shared_ptr<Field> lon(new Field("lon", {dim2}, {"ncol"}, nparts, 0));
  const Real lat_data[] = {-0.5, 0.5, 0.5, -0.5};
  const Real lon_data[] = {-0.5, -0.5, 0.5, 0.5};
  for (int i = 0; i < nparts; ++i) {
    lat->set_part_size(i, part_size);
    lat->set_part_data(i, &lat_data[part_size*i]);
    lon->set_part_size(i, part_size);
    lon->set_part_data(i, &lon_data[part_size*i]);
  }
  lat->commit();
  lon->commit();

  // Allocate dummy lat/lon
  const int dum_nparts = nparts-1;
  const int dum_part_size = dim2 / dum_nparts;
  std::shared_ptr<Field> dum_lat(new Field("dum_lat", {dim2}, {"ncol"}, dum_nparts, 0));
  std::shared_ptr<Field> dum_lon(new Field("lon", {dim2}, {"ncol"}, dum_nparts, 0));
  for (int i = 0; i < dum_nparts; ++i) {
    dum_lat->set_part_size(i, dum_part_size);
    dum_lat->set_part_data(i, &lat_data[dum_part_size*i]);
    dum_lon->set_part_size(i, dum_part_size);
    dum_lon->set_part_data(i, &lon_data[dum_part_size*i]);
  }
  dum_lat->commit();
  dum_lon->commit();

  // Test bounding box field
  // bounds set to lat = {0.0, 1.0}, lon = {0.0, 1.0}
  const auto bounding_box_sname = "bounding_box";
  auto bb_stat = create_stat(bounding_box_sname, comm);
  auto bounding_box_stat = dynamic_cast<FieldBoundingBox *>(bb_stat.get());
  REQUIRE_THROWS(bounding_box_stat->compute(foo)); // initialize() required
  REQUIRE_THROWS(bounding_box_stat->initialize(dum_lat, dum_lon)); // dum_lat wrong name
  bounding_box_stat->initialize(lat, dum_lon);
  REQUIRE_THROWS(bounding_box_stat->compute(foo)); // dum_lon wrong size
  bounding_box_stat->initialize(lat, lon);
  const auto bounding_box_field = bounding_box_stat->compute(foo);
  const Real bounding_box_expected[] = {
      0.0, 0.0, 12.0, 0.0,
      0.0, 0.0, 14.0, 0.0,
      0.0, 0.0, 16.0, 0.0,
      0.0, 0.0, 18.0, 0.0,
      0.0, 0.0, 20.0, 0.0,
      0.0, 0.0, 22.0, 0.0};
  for (int i = 0; i < dim0*dim1*dim2; ++i)
    REQUIRE (bounding_box_expected[i]==bounding_box_field.data<Real>()[i]);

  // Allocate area
  std::shared_ptr<Field> area(new Field("area", {dim2}, {"ncol"}, nparts, 0));
  const Real area_data[] = {0.5, 1.0, 1.5, 2.0};
  for (int i = 0; i < nparts; ++i) {
    area->set_part_size(i, part_size);
    area->set_part_data(i, &area_data[part_size*i]);
  }
  area->commit();

  // Allocate dummy area
  std::shared_ptr<Field> dum_area(new Field("area", {dim2}, {"ncol"}, dum_nparts, 0));
  for (int i = 0; i < dum_nparts; ++i) {
    dum_area->set_part_size(i, dum_part_size);
    dum_area->set_part_data(i, &area_data[dum_part_size*i]);
  }
  dum_area->commit();

  // Test zonal mean
  // bounds set to lat = {0.0, 1.0}
  const auto zonal_mean_sname = "zonal_mean";
  auto zm_stat = create_stat(zonal_mean_sname, comm);
  auto zonal_mean_stat = dynamic_cast<FieldZonalMean *>(zm_stat.get());
  REQUIRE_THROWS(zonal_mean_stat->compute(foo)); // initialize() required
  REQUIRE_THROWS(zonal_mean_stat->initialize(dum_lat, dum_area)); // dum_lat wrong name
  REQUIRE_THROWS(zonal_mean_stat->initialize(lat, dum_area)); // dum_area wrong size
  zonal_mean_stat->initialize(lat, area);
  const auto zonal_mean_field = zonal_mean_stat->compute(foo);
  const Real zonal_area = area_data[1] + area_data[2];
  const Real zonal_mean_expected[] = {
    ((1.0 * area_data[1]) + (12.0 * area_data[2])) / zonal_area,
    ((3.0 * area_data[1]) + (14.0 * area_data[2])) / zonal_area,
    ((5.0 * area_data[1]) + (16.0 * area_data[2])) / zonal_area,
    ((7.0 * area_data[1]) + (18.0 * area_data[2])) / zonal_area,
    ((9.0 * area_data[1]) + (20.0 * area_data[2])) / zonal_area,
    ((11.0 * area_data[1]) + (22.0 * area_data[2])) / zonal_area};
  for (int i = 0; i < dim0*dim1; ++i)
    REQUIRE (zonal_mean_expected[i]==zonal_mean_field.data<Real>()[i]);
}
