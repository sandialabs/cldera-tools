#include "profiling/control/cldera_control.hpp"

#include <catch2/catch.hpp>

TEST_CASE ("control") {
  using namespace cldera;

  SECTION ("basic") {
    REQUIRE (true);
  }

  // Check uniform alloc size 
  SECTION ("part_alloc_size") {
    REQUIRE (true);
  }

  SECTION ("update") {
    //REQUIRE_THROWS(yr.update(b1,1,1)); // layout mismatch
    REQUIRE (true);
  }
}
