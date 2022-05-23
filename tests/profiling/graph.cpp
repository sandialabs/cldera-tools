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

  // A graph from an input file
  {
    std::string filename = "./cldera_graph_input.yaml";
    
    cldera::GraphFactory factory(filename, true);

    std::shared_ptr<cldera::Graph> graph = factory.build_graph();

    REQUIRE(graph->get_num_vertices() == 4);
    REQUIRE(graph->get_num_edges() == 3);
    REQUIRE(graph->is_cyclic() == false);
  }

  // A graph with a cycle A->B, B->A
  {
    using cldera::GraphVertex;
    std::map<std::string, std::shared_ptr<GraphVertex> > graph_vertices 
      = { {"A",std::make_shared<GraphVertex>(GraphVertex("A"))}, {"B",std::make_shared<GraphVertex>(GraphVertex("B"))} };
    std::map<std::string, std::vector<std::string> > graph_edges = { {"A",{"B"}}, {"B",{"A"}} };
    cldera::Graph graph(graph_vertices, graph_edges);

    REQUIRE(graph.get_num_vertices() == 2);
    REQUIRE(graph.get_num_edges() == 2);
    REQUIRE(graph.is_cyclic() == true);
    REQUIRE(graph.get_children("A")[0] == "B");
    REQUIRE(graph.get_parents("A")[0] == "B");
    REQUIRE(graph.get_children("B")[0] == "A");
    REQUIRE(graph.get_parents("B")[0] == "A");
  }

  // A graph with a cycle A->B, B->C, C->A,D
  {
    using cldera::GraphVertex;
    std::map<std::string, std::shared_ptr<GraphVertex> > graph_vertices 
      = { {"A",std::make_shared<GraphVertex>(GraphVertex("A"))}, {"B",std::make_shared<GraphVertex>(GraphVertex("B"))}, 
          {"C",std::make_shared<GraphVertex>(GraphVertex("C"))}, {"D",std::make_shared<GraphVertex>(GraphVertex("D"))} };
    std::map<std::string, std::vector<std::string> > graph_edges = { {"A",{"B"}}, {"B",{"C"}}, {"C",{"D","A"}} };
    cldera::Graph graph(graph_vertices, graph_edges);

    REQUIRE(graph.get_num_vertices() == 4);
    REQUIRE(graph.get_num_edges() == 4);
    REQUIRE(graph.is_cyclic() == true);
    REQUIRE(graph.get_children("A")[0] == "B");
    REQUIRE(graph.get_children("B")[0] == "C");
    REQUIRE(graph.get_parents("C")[0] == "B");
    REQUIRE(graph.get_parents("B")[0] == "A");
    REQUIRE(graph.get_children_recursive("A",-1).size() == 4);
    REQUIRE(graph.get_children_recursive("A",1).size() == 1);
  }

}
