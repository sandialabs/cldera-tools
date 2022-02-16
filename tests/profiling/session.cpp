#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_session.hpp"
#include "cldera_config.h"

TEST_CASE ("session") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto& s = ProfilingSession::instance();

#ifdef CLDERA_DEBUG
  // Not yet inited, so can't clean up
  REQUIRE_THROWS (s.clean_up());
#endif

  s.init (comm);

  ekat::any foo,bar;

  // Create a simple 'any', set it in the session
  foo.reset<double>(1.0);
  s.set_data("foo",foo);

  // Cannot overwrite.
  REQUIRE_THROWS (s.set_data("foo",foo));

  // Check data is in there, has correct type, and correct value
  REQUIRE_NOTHROW( bar = s.get_data("foo") );
  REQUIRE_NOTHROW (ekat::any_cast<double>(bar));
  REQUIRE (ekat::any_cast<double>(bar)==1.0);

#ifdef CLDERA_DEBUG
  // Cannot double init
  REQUIRE_THROWS(s.init (comm));
#endif

  s.clean_up();
}
