#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_archive.hpp"
#include "profiling/stats/cldera_field_global_max.hpp"

TEST_CASE ("archive") {
  using namespace cldera;

  const ekat::Comm comm(MPI_COMM_WORLD);

  int ymd = 20220915;
  int tod = 43000;
  TimeStamp ts(ymd,tod);

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix","archive_tests");
  params.set("Flush Frequency",5);

  ProfilingArchive archive(comm,ts,ts,params);

  std::vector<Real> foo_data (20,3.0);
  archive.add_field(Field("foo",{5,4},{"col","lev"},foo_data.data()));
  REQUIRE (archive.has_field("foo"));
  REQUIRE_THROWS (archive.get_field("bar"));

  auto foo = archive.get_field("foo");
  REQUIRE (foo.layout().size()==20);

  FieldGlobalMax foo_max(comm,ekat::ParameterList("foo_global_max"));
  foo_max.set_field(foo);
  foo_max.create_stat_field();
  archive.update_stat("foo",foo_max.name(),foo_max.compute(ts));
  archive.end_timestep(ts+=86400);
}
