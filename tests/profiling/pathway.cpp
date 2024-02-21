#include <catch2/catch.hpp>

#include "profiling/stats/cldera_register_stats.hpp"
#include "profiling/cldera_profiling_archive.hpp"
#include "profiling/cldera_graph_factory.hpp"
#include "profiling/cldera_pathway_factory.hpp"
#include "profiling/cldera_profiling_types.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <memory>
#include <vector>

TEST_CASE ("pathway") {

  ekat::Comm comm(MPI_COMM_WORLD);
  cldera::TimeStamp time = {19910701, 3600};

  cldera::register_stats();

  // A Pathway via PathwayFactory
  {
    // initialize simple field
    constexpr int foo_size = 5;
    const std::vector<int> foo_sizes(1, foo_size);
    std::vector<cldera::Real> foo_data(foo_size);
    std::vector<std::string> dimnames = {"mydim"};
    cldera::Field foo("foo", foo_sizes, dimnames, foo_data.data());

    // initialize another simple field, the same size as foo
    std::vector<cldera::Real> bar_data(foo_size);
    cldera::Field bar("bar", foo_sizes, dimnames, bar_data.data());

    // Create an archive
    cldera::ProfilingArchive archive (comm,time,time,{});
    archive.add_field(foo);
    archive.add_field(bar);

    // load the input deck and construct the pathway
    auto params = ekat::parse_yaml_file("./cldera_pathway_input.yaml");
    cldera::PathwayFactory pathway_factory(params, comm, true);
    auto pathway = pathway_factory.build_pathway(archive, std::cout);

    // check that all pathway tests pass
    std::iota(foo_data.begin(), foo_data.end(), 1.0);
    std::iota(bar_data.begin(), bar_data.end(), 1.0);
    REQUIRE(pathway->run_pathway_tests(comm, time));

    foo_data[2] = 8.0;
    time = {19910701, 7200};
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    foo_data[2] = -1.0;
    time = {19910701, 10800};
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    foo_data[2] = 1.0;
    bar_data[2] = 2.0;
    time = {19910701, 14400};
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    bar_data[2] = -1.0;
    time = {19910701, 18000};
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    bar_data[2] = 1.0;
    time = {19910701, 21600};
    //REQUIRE(pathway->run_pathway_tests(comm, time)); // TODO: is this a bug?

    pathway->dump_test_history_to_yaml("test_history.yaml");
  }

}
