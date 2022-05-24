#ifndef CLDERA_GRAPH_HPP
#define CLDERA_GRAPH_HPP

#include "cldera_graph_vertex.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <memory>  // for std::shared_ptr
#include <map>     // for std::map

namespace cldera
{

/**
 * A class that stores simple information such as directed edges and vertices.
 * It is hardcoded to a GraphVertex class, but should be easy to support 
 * templates or vertex polymorphism later.
 */
class Graph {
public:
  Graph (std::map<std::string, std::shared_ptr<GraphVertex> > vertices,
         std::map<std::string, std::vector<std::string> > edges)
  : m_vertices(vertices),
    m_edges(edges)
  {}

  // Print a dot graph structure to stdout which can then be used to generate an image
  void generate_dot_graph(std::ostream& out = std::cout) const;

  // Computes whether the graph is cyclic or not
  bool is_cyclic();

  // Setters

  // None at the moment


  // Getters

  // Get the vertex object keyed by the string name
  std::shared_ptr<GraphVertex> get_vertex(const std::string& name) { return m_vertices[name]; }

  // Get the names of children of a node keyed by the string name
  std::vector<std::string> get_children(const std::string& name) { return m_edges[name]; }

  // Get the names of children of a node keyed by the string name by traversing a given depth
  std::vector<std::string> get_children_recursive(const std::string& name, const int recursion_depth = -1);

  // Get the names of parents of a given node keyed by the string name
  std::vector<std::string> get_parents(const std::string& name);

  // Return the size of m_vertices
  size_t get_num_vertices() const { return m_vertices.size(); }

  // Return the number of edges, which is the sum of m_edges[key].size() for each valid key
  size_t get_num_edges() const;

private:
  // The vertices of the graph, which is a map that takes names and returns GraphVertex objects
  std::map<std::string, std::shared_ptr<GraphVertex> > m_vertices;

  // The *directed* edges of the graph. A *directed* edge is the pair of vertices defined by m_edges[SOURCE]=DESTINATION
  std::map<std::string, std::vector<std::string> > m_edges;
};

} // namespace cldera

#endif // CLDERA_GRAPH_HPP
