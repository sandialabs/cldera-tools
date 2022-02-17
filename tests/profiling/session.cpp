#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_session.hpp"
#include "cldera_config.h"

#include "profiling/cldera_profiling_types.hpp"

TEST_CASE ("session") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto& s = ProfilingSession::instance();

#ifdef CLDERA_DEBUG
  // Not yet inited, so can't clean up
  REQUIRE_THROWS (s.clean_up());
#endif

  s.init (comm);

  // Create a simple 'any', set it in the session
  s.create<double>("foo",1.0);

  // Cannot overwrite.
  REQUIRE_THROWS (s.create<double>("foo",1.0));

  // But you can use create_or_get
  REQUIRE_NOTHROW (s.create_or_get<double>("foo"));

  // Cannot get/remove what's not there
  REQUIRE_THROWS (s.get<int>("bar"));
  REQUIRE_THROWS (s.remove("bar"));

  // Check data is stored, has correct type, and correct value
  REQUIRE_THROWS ( s.get<int>("foo") );
  REQUIRE (s.get<double>("foo")==1.0);

#ifdef CLDERA_DEBUG
  // Cannot double init
  REQUIRE_THROWS(s.init (comm));
#endif

  s.clean_up();
}
