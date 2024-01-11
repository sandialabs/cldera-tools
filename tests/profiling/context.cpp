#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_context.hpp"
#include "cldera_config.h"

#include "profiling/cldera_profiling_types.hpp"

TEST_CASE ("session") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  ProfilingContext c("test");

#ifdef CLDERA_DEBUG
  // Not yet inited, so can't clean up
  REQUIRE_THROWS (c.clean_up());
#endif

  c.init (comm,{});

#ifdef CLDERA_DEBUG
  // Cannot double init
  REQUIRE_THROWS(c.init (comm,{}));
#endif

  // Create a simple 'any', set it in the session
  c.create<double>("foo",1.0);

  // Cannot overwrite.
  REQUIRE_THROWS (c.create<double>("foo",1.0));

  // Cannot get what's not there
  REQUIRE_THROWS (c.get<int>("bar"));

  // Check data is stored, has correct type, and correct value
  REQUIRE_THROWS ( c.get<int>("foo") );
  REQUIRE (c.get<double>("foo")==1.0);

  c.clean_up();
}
