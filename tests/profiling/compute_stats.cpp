#include "profiling/stats/cldera_register_stats.hpp"
#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <catch2/catch.hpp>

#include <map>

TEST_CASE ("stats_tests") {
  using namespace cldera;

  // REQUIRE_THROWS causes some start_timer calls to not be matched
  // by a corresponding stop_timer. When that happens, a subsequent
  // call to start_timer for that timer would throw. For this test,
  // we can just disable timings
  timing::TimingSession::instance().toggle_session(false);

  // Register all stats in the factory
  register_stats();

  ekat::Comm comm(MPI_COMM_WORLD);

  int ymd = 20220915;
  int tod = 43000;
  TimeStamp time(ymd,tod);

  SECTION ("global_stats") {
    constexpr int dim0 = 2;
    constexpr int dim1 = 3;

    // Create raw data
    double raw_data[dim0*dim1];
    std::iota(raw_data,raw_data+dim0*dim1,0);

    Field f("f",{dim0,dim1},{"dim0","dim1"});
    f.set_data(raw_data);

    // Compute all supported stats
    auto stats_types = {"global_max", "global_min", "global_sum", "global_avg"};
    std::map<std::string,Field> stat_fields;
    std::map<std::string,Field> expected;
    for (std::string type : stats_types) {
      auto sname = "my_" + type;
      ekat::ParameterList pl (sname);
      pl.set("Type",type);
      auto stat = StatFactory::instance().create(type,comm,pl);
      REQUIRE_THROWS (stat->compute(time)); // field not set yet
      stat->set_field(f);
      stat->create_stat_field();
      stat_fields.emplace(sname,stat->compute(time));

      expected.emplace(sname,Field(sname,stat_fields.at(sname).layout(),DataAccess::Copy));
      expected.at(sname).commit();
    }

    auto stats_names = {"my_global_max", "my_global_min", "my_global_sum", "my_global_avg"};
    expected["my_global_max"].data_nonconst<Real>()[0] = 5;
    expected["my_global_min"].data_nonconst<Real>()[0] = 0;
    expected["my_global_sum"].data_nonconst<Real>()[0] = 15;
    expected["my_global_avg"].data_nonconst<Real>()[0] = 15.0/6;

    for (auto sname : stats_names) {
      REQUIRE (expected.at(sname).data<Real>()[0]==stat_fields.at(sname).data<Real>()[0]);
    }
  }

  SECTION ("stat_int_field") {
    constexpr int dim0 = 10;

    // Create raw data
    int raw_data[dim0];
    std::iota(raw_data,raw_data+dim0,1);

    Field f("f",{dim0},{"dim0"},DataAccess::View,IntType);
    f.set_data(raw_data);

    auto stat = StatFactory::instance().create("global_sum",comm,ekat::ParameterList("global_sum"));
    stat->set_field(f);
    stat->create_stat_field();
    const auto computed = stat->compute(time).data<int>()[0];
    const auto expected = (dim0+1)*dim0/2;
    REQUIRE (computed==expected);
  }

  SECTION ("stats along columns") {
    constexpr int dim0 = 2;
    constexpr int dim1 = 3;
    constexpr int dim2 = 4;

    // Create raw data
    double raw_data[dim0*dim1*dim2];
    std::iota(raw_data,raw_data+dim0*dim1*dim2,0);

    Field f("f", {dim0,dim1,dim2}, {"lev", "dim", "ncol"});
    f.set_data(raw_data);

    auto stats_names = {"max_along_columns", "min_along_columns",
        "sum_along_columns", "avg_along_columns"};
    std::map<std::string, Field> stat_fields;
    std::map<std::string, Field> expected;
    for (auto sname : stats_names) {
      auto stat = StatFactory::instance().create(sname,comm,ekat::ParameterList(sname));
      stat->set_field(f);
      stat->create_stat_field();
      stat_fields.emplace(sname, stat->compute(time));
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

  SECTION ("stats along columns with parts") {
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
      foo.set_part_extent(i, part_size);
      foo.set_part_data(i, foo_data[i].data());
    }
    foo.commit();

    // Compute all supported stats along columns
    auto stats_names = {"max_along_columns", "min_along_columns",
        "sum_along_columns", "avg_along_columns"};
    std::map<std::string, Field> stat_fields;
    std::map<std::string, Field> expected;
    for (auto sname : stats_names) {
      auto stat = StatFactory::instance().create(sname,comm,ekat::ParameterList(sname));
      stat->set_field(foo);
      stat->create_stat_field();
      stat_fields.emplace(sname, stat->compute(time));
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

  SECTION ("stats_with_bounds") {
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
      foo.set_part_extent(i, part_size);
      foo.set_part_data(i, foo_data[i].data());
    }
    foo.commit();

    // Allocate lat/lon
    Field lat("lat", {dim2}, {"ncol"}, nparts, 0);
    Field lon("lon", {dim2}, {"ncol"}, nparts, 0);
    const Real lat_data[] = {-0.5, 0.5, 0.5, -0.5};
    const Real lon_data[] = {-0.5, -0.5, 0.5, 0.5};
    for (int i = 0; i < nparts; ++i) {
      lat.set_part_extent(i, part_size);
      lat.set_part_data(i, &lat_data[part_size*i]);
      lon.set_part_extent(i, part_size);
      lon.set_part_data(i, &lon_data[part_size*i]);
    }
    lat.commit();
    lon.commit();

    // Allocate dummy lat/lon
    const int dum_nparts = nparts-1;
    const int dum_part_size = dim2 / dum_nparts;
    Field dum_lat("dum_lat", {dim2}, {"ncol"}, dum_nparts, 0);
    Field dum_lon("lon", {dim2}, {"ncol"}, dum_nparts, 0);
    for (int i = 0; i < dum_nparts; ++i) {
      dum_lat.set_part_extent(i, dum_part_size);
      dum_lat.set_part_data(i, &lat_data[dum_part_size*i]);
      dum_lon.set_part_extent(i, dum_part_size);
      dum_lon.set_part_data(i, &lon_data[dum_part_size*i]);
    }
    dum_lat.commit();
    dum_lon.commit();

    // Test bounded field
    SECTION("bounded") {
      auto bounded_pl = ekat::ParameterList("bounded");
      bounded_pl.set<std::vector<Real>>("Bounds", {12.1, 20.1});
      const auto bounded_stat = StatFactory::instance().create("bounded",comm,bounded_pl);
      bounded_stat->set_field(foo);
      bounded_stat->create_stat_field();
      const auto bounded_field = bounded_stat->compute(time);
      const Real bounded_expected[] = {
          0.0, 0.0, 0.0, 13.0,
          0.0, 0.0, 14.0, 15.0,
          0.0, 0.0, 16.0, 17.0,
          0.0, 0.0, 18.0, 19.0,
          0.0, 0.0, 20.0, 0.0,
          0.0, 0.0, 0.0, 0.0};
      for (int i = 0; i < dim0*dim1*dim2; ++i)
        REQUIRE (bounded_expected[i]==bounded_field.data<Real>()[i]);
    }

    // Test bounding box field
    SECTION ("bounding_box") {
      auto bounding_box_pl = ekat::ParameterList("bounding_box");
      bounding_box_pl.set<std::vector<Real>>("Latitude Bounds", {0.0, 1.0});
      bounding_box_pl.set<std::vector<Real>>("Longitude Bounds", {0.0, 1.0});
      auto bounding_box_stat = StatFactory::instance().create("bounding_box",comm,bounding_box_pl);
      bounding_box_stat->set_field(foo);
      REQUIRE_THROWS(bounding_box_stat->create_stat_field()); // set_aux_fields() required
      REQUIRE_THROWS(bounding_box_stat->set_aux_fields(dum_lat, dum_lon)); // dum_lat wrong name
      REQUIRE_THROWS(bounding_box_stat->set_aux_fields(lat, dum_lon));  // dum_lon wrong nparts
      bounding_box_stat->set_aux_fields(lat, lon);
      bounding_box_stat->create_stat_field();
      const auto bounding_box_field = bounding_box_stat->compute(time);
      const Real bounding_box_expected[] = {
          0.0, 0.0, 12.0, 0.0,
          0.0, 0.0, 14.0, 0.0,
          0.0, 0.0, 16.0, 0.0,
          0.0, 0.0, 18.0, 0.0,
          0.0, 0.0, 20.0, 0.0,
          0.0, 0.0, 22.0, 0.0};
      for (int i = 0; i < dim0*dim1*dim2; ++i)
        REQUIRE (bounding_box_expected[i]==bounding_box_field.data<Real>()[i]);
    }

    // Test zonal mean
    SECTION ("zonal_mean") {
      // Allocate area
      Field area("area", {dim2}, {"ncol"}, nparts, 0);
      const Real area_data[] = {0.5, 1.0, 1.5, 2.0};
      for (int i = 0; i < nparts; ++i) {
        area.set_part_extent(i, part_size);
        area.set_part_data(i, &area_data[part_size*i]);
      }
      area.commit();

      // Allocate dummy area
      Field dum_area("area", {dim2}, {"ncol"}, dum_nparts, 0);
      for (int i = 0; i < dum_nparts; ++i) {
        dum_area.set_part_extent(i, dum_part_size);
        dum_area.set_part_data(i, &area_data[dum_part_size*i]);
      }
      dum_area.commit();

      auto zonal_mean_pl = ekat::ParameterList("zonal_mean");
      zonal_mean_pl.set<std::vector<Real>>("Latitude Bounds", {0.0, 1.0});
      auto zonal_mean_stat = StatFactory::instance().create("zonal_mean",comm,zonal_mean_pl);
      zonal_mean_stat->set_field(foo);
      REQUIRE_THROWS(zonal_mean_stat->create_stat_field()); // set_aux_fields() required
      REQUIRE_THROWS(zonal_mean_stat->set_aux_fields(dum_lat, dum_area)); // dum_lat wrong name
      REQUIRE_THROWS(zonal_mean_stat->set_aux_fields(lat, dum_area)); // dum_area wrong size
      zonal_mean_stat->set_aux_fields(lat, area);
      zonal_mean_stat->create_stat_field();
      const auto zonal_mean_field = zonal_mean_stat->compute(time);
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
  }
}
