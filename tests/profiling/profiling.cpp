#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_session.hpp"
#include "cldera_config.h"

TEST_CASE ("dummy") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto& s = ProfilingSession::instance();

#ifdef CLDERA_DEBUG
  // Not yet inited, so can't clean up
  REQUIRE_THROWS (s.clean_up());
#endif

  s.init (comm);

#ifdef CLDERA_DEBUG
  // Cannot double init
  REQUIRE_THROWS(s.init (comm));
#endif

  s.clean_up();
}
