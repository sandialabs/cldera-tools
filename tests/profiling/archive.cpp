#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_archive.hpp"

TEST_CASE ("archive") {
  using namespace cldera;

  ProfilingArchive archive;

  std::vector<Real> foo_data (20);
  archive.add_field(Field("foo",{5,4},{"col","lev"},foo_data.data()));

  REQUIRE (archive.has_field("foo"));
  REQUIRE (archive.get_field("foo").layout().size()==20);
  REQUIRE_THROWS (archive.get_field("bar"));
  REQUIRE_THROWS (archive.get_stat_history("bar",StatType::Max));
  REQUIRE (archive.get_stat_history("foo",StatType::Max).size()==0);
}
