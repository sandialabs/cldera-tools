#include <catch2/catch.hpp>

#include "profiling/cldera_bounds_field_test.hpp"
#include "profiling/cldera_field.hpp"
#include "profiling/cldera_profiling_test_manager.hpp"
#include "profiling/cldera_profiling_types.hpp"


TEST_CASE ("test_manager") {
  using namespace cldera;

  // Initialize test manager
  ProfilingTestManager test_manager;

  // Initialize simple field
  constexpr int foo_size = 5;
  const std::vector<int> foo_sizes(1, foo_size);
  std::vector<Real> foo_data(foo_size);
  const auto foo = std::make_shared<const Field>("foo", foo_sizes, foo_data.data());

  // Initialize simple field test
  const std::string field_test_name = "Test bounds of foo";
  const Real min = 0.0, max = 6.0;
  const Bounds bounds{min, max};
  const auto field_test = std::make_shared<BoundsFieldTest>(field_test_name, foo, bounds);

  // Test add field test to test manager
  REQUIRE_THROWS(test_manager.add_field_test(nullptr)); // Throw on nullptr
  REQUIRE_NOTHROW(test_manager.add_field_test(field_test));
  REQUIRE_THROWS(test_manager.add_field_test(field_test)); // Throw if field_test already exists
  
  // Test has_field_test
  REQUIRE(test_manager.has_field_test("Test does not exist") == false);
  REQUIRE(test_manager.has_field_test(field_test_name) == true);

  // Test field test runner
  const ekat::Comm comm(MPI_COMM_WORLD);
  REQUIRE_THROWS(test_manager.run_field_test("Test does not exist", comm));
  std::iota(foo_data.begin(), foo_data.end(), 1); // field within bounds
  REQUIRE(test_manager.run_field_test(field_test_name, comm) == true);
  foo_data[2] = -1.5; // field min out of bounds
  REQUIRE(test_manager.run_field_test(field_test_name, comm) == false);
}
