#include <catch2/catch.hpp>

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

  // A Pathway via PathwayFactory
  {
    // initialize simple field
    constexpr int foo_size = 5;
    const std::vector<int> foo_sizes(1, foo_size);
    std::vector<cldera::Real> foo_data(foo_size);
    std::vector<std::string> dimnames = {"mydim"};
    const auto foo = std::make_shared<const cldera::Field>("foo", foo_sizes, dimnames, foo_data.data());

    // initialize another simple field, the same size as foo
    std::vector<cldera::Real> bar_data(foo_size);
    const auto bar = std::make_shared<const cldera::Field>("bar", foo_sizes, dimnames, bar_data.data());

    std::map<std::string, std::shared_ptr<const cldera::Field> > fields;
    fields["foo"] = foo;
    fields["bar"] = bar;

    // load the input deck and construct the pathway
    std::string filename = "./cldera_pathway_input.yaml";
    cldera::PathwayFactory pathway_factory(filename, fields, true);
    std::shared_ptr<cldera::Pathway> pathway = pathway_factory.build_pathway(std::cout);
    cldera::TimeStamp time = {1991, 3600};

    // check that all pathway tests pass
    std::iota(foo_data.begin(), foo_data.end(), 1.0);
    std::iota(bar_data.begin(), bar_data.end(), 1.0);
    REQUIRE(pathway->run_pathway_tests(comm, time));

    foo_data[2] = 8.0;
    time.tod = 7200;
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    foo_data[2] = -1.0;
    time.tod = 10800;
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    foo_data[2] = 1.0;
    bar_data[2] = 2.0;
    time.tod = 14400;
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    bar_data[2] = -1.0;
    time.tod = 18000;
    REQUIRE(!pathway->run_pathway_tests(comm, time));

    bar_data[2] = 1.0;
    time.tod = 21600;
    //REQUIRE(pathway->run_pathway_tests(comm, time)); // TODO: is this a bug?

    pathway->dump_test_history_to_yaml("test_history.yaml");
  }

}
