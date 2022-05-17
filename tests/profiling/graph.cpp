#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_archive.hpp"
#include "profiling/cldera_graph_factory.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <memory>
#include <vector>

TEST_CASE ("graph") {

  std::string filename = "./cldera_graph_input.yaml";
  
  cldera::GraphFactory factory(filename, true);

  std::shared_ptr<cldera::Graph> graph = factory.build_graph();

  REQUIRE(graph->get_num_vertices() == 4);
  REQUIRE(graph->get_num_edges() == 3);

}
