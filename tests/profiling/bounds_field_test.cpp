#include <catch2/catch.hpp>

#include "profiling/cldera_bounds_field_test.hpp"
#include "profiling/cldera_field.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include "memory"

TEST_CASE ("bounds_field_test") {
  using namespace cldera;

  // Initialize simple field
  constexpr int foo_size = 5;
  const std::vector<int> foo_sizes(1, foo_size);
  std::vector<Real> foo_data(foo_size);
  std::vector<std::string> dimnames = {"mydim"};
  const auto foo = std::make_shared<const Field>("foo", foo_sizes, dimnames, foo_data.data());

  // Initialize bounds field test
  const std::string bounds_field_test_name = "Test bounds of foo";
  const Real min = 0.0, max = 6.0;
  const Bounds bounds{min, max};
  const BoundsFieldTest bounds_field_test(bounds_field_test_name, foo, bounds);

  // Check name
  REQUIRE(bounds_field_test.name() == bounds_field_test_name);

  // Check if field is within bounds
  std::iota(foo_data.begin(), foo_data.end(), 1);
  const ekat::Comm comm(MPI_COMM_WORLD);
  REQUIRE(bounds_field_test.test(comm) == true);

  // Check min failure
  foo_data[2] = -1.5;
  REQUIRE(bounds_field_test.test(comm) == false);

  // Check max failure
  foo_data[2] = 7.5;
  REQUIRE(bounds_field_test.test(comm) == false);
}
