#include "timing/cldera_timing_session.hpp"

#include <ekat/ekat_assert.hpp>

#include <catch2/catch.hpp>
#include <iostream>

TEST_CASE ("timing")
{
  using namespace cldera::timing;

  ekat::Comm comm(MPI_COMM_WORLD);
  const int rank = comm.rank();
  const int size = comm.size();

  auto& s = TimingSession::instance();

  s.start_timer("first");
  REQUIRE_THROWS(s.start_timer("first"));
  s.stop_timer("first");
  REQUIRE_THROWS(s.stop_timer("first"));
  s.start_timer("first");
  s.start_timer("second");
  s.stop_timer("second");
  s.stop_timer("first");

  s.dump(std::cout,comm);
}
