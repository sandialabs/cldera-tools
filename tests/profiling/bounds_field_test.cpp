#include <catch2/catch.hpp>

#include "profiling/cldera_bounds_field_test.hpp"
#include "profiling/cldera_field_test_factory.hpp"
#include "profiling/cldera_field.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <map>
#include "memory"

TEST_CASE ("bounds_field_test") {
  using namespace cldera;

  const ekat::Comm comm(MPI_COMM_WORLD);

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
  const BoundsFieldTest bounds_field_test(bounds_field_test_name, foo, bounds, comm);
  TimeStamp time = {1,1};

  // Test BoundsFieldTest
  {
    BoundsFieldTest bounds_field_test(bounds_field_test_name, foo, bounds, comm);
    bounds_field_test.set_save_history(true);
    
    // Check name
    REQUIRE(bounds_field_test.name() == bounds_field_test_name);

    // Check if field is within bounds
    std::iota(foo_data.begin(), foo_data.end(), 1.0);
    REQUIRE(bounds_field_test.test(time) == true);

    // Check min failure
    foo_data[2] = -1.5;
    REQUIRE(bounds_field_test.test(time) == false);

    // Check max failure
    foo_data[2] = 7.5;
    REQUIRE(bounds_field_test.test(time) == false);
    
    // Check if history has been saved for both failures
    REQUIRE(bounds_field_test.get_test_history().size() == 2);
  }

  // Test BoundsFieldTest, MinFieldTest, and MaxFieldTest via FieldTestFactory
  {
    // Create a map of fields
    std::map<std::string, std::shared_ptr<const Field> > fields;
    fields.emplace("foo",foo);

    // Create a field test factory
    FieldTestFactory field_test_factory("cldera_field_test_input.yaml", fields, comm, true);
    std::map<std::string, std::shared_ptr<FieldTest> > tests = field_test_factory.build_field_tests(std::cout);
    
    // Check that the name used to access the test is the same as its member
    REQUIRE(tests["foobounds"]->name() == "foobounds");
    REQUIRE(tests["foomax"]->name() == "foomax");
    REQUIRE(tests["foomin"]->name() == "foomin");

    // Check if field is within bounds
    std::iota(foo_data.begin(), foo_data.end(), 1.0);
    REQUIRE(tests["foobounds"]->test(time) == true);
    REQUIRE(tests["foomax"]->test(time) == true);
    REQUIRE(tests["foomin"]->test(time) == true);

    // Check min failure
    foo_data[2] = -1.5;
    REQUIRE(tests["foobounds"]->test(time) == false);
    REQUIRE(tests["foomax"]->test(time) == true);
    REQUIRE(tests["foomin"]->test(time) == false);

    // Check max failure
    foo_data[2] = 7.5;
    REQUIRE(tests["foobounds"]->test(time) == false);
    REQUIRE(tests["foomax"]->test(time) == false);
    REQUIRE(tests["foomin"]->test(time) == true);
  }
}
