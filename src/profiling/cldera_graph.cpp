#include "cldera_graph.hpp"

namespace cldera
{

size_t Graph::get_num_edges() const
{ 
  size_t count = 0;
  for (auto& pair : m_edges)
    count += pair.second.size();
  
  return count;
}

void Graph::generate_dot_graph(std::ostream& out) const
{
  out << "digraph graphname {\n";

  for(auto& pair : m_edges)
  {
    for (const auto& dst : pair.second)
    {
      out << pair.first << " -> " << dst << "\n";
    }
  }
  out << "}\n";
}

/**
 * Find all children of a given vertex named name. 
 * recursion_depth = -1 is infinite by default, 
 * recursion_depth = 0 is empty, and 
 * recursion_depth = 1 is just the children of name
 */
std::vector<std::string> Graph::get_children_recursive(const std::string& name, const int recursion_depth)
{
  EKAT_REQUIRE_MSG (recursion_depth >= -1, "Error! Graph.get_children_recursive: recursion_depth must be -1 or greater!\n");
  std::vector<std::string> children;

  int depth = 0;

  std::vector<std::string> vertices_to_traverse = { name };
  // As long as there are children to traverse and the depth has not been exceeded
  while ((recursion_depth == -1 || depth < recursion_depth) && vertices_to_traverse.size() > 0)
  {
    // grab the children of each node to traverse and insert them to the list
    std::vector<std::string> children_to_traverse;
    for (size_t i = 0; i < vertices_to_traverse.size(); ++i)
    {
      std::vector<std::string> vertex_children = m_edges[vertices_to_traverse[i]];
      children_to_traverse.insert(children_to_traverse.end(), std::make_move_iterator(vertex_children.begin()), std::make_move_iterator(vertex_children.end()));
    }

    // Check for duplicate indices (i.e. cycles)
    std::vector<size_t> duplicate_indices;
    for (size_t i = 0; i < children_to_traverse.size(); ++i)
    {
      for (size_t j = 0; j < children.size(); ++j)
      {
        if (children_to_traverse[i] == children[j])
        {
          duplicate_indices.push_back(i);
        }
      }
    }

    // Remove duplicates in reverse order since it resizes the array
    for (size_t i = 0; i < duplicate_indices.size(); ++i)
    {
      children_to_traverse.erase(children_to_traverse.begin() + duplicate_indices[duplicate_indices.size() - 1 - i]);
    }

    // Add all non-duplicates to the list of children
    children.insert(children.end(), children_to_traverse.begin(), children_to_traverse.end());

    // Prepare for the next depth search
    vertices_to_traverse = children_to_traverse;
    depth++;
  }
  
  return children;
}

std::vector<std::string> Graph::get_parents(const std::string& name)
{
  std::vector<std::string> parents;

  // Loop over potential parent vertices
  for (auto& pair : m_vertices)
  {
    std::string parent_name = pair.first;
    std::vector<std::string> children = m_edges[parent_name];

    // Loop over children of potential parents
    for(size_t i = 0; i < children.size(); ++i)
    {
      // If the potential parent has a child named name, it's a parent of name
      if(children[i] == name)
      {
        parents.push_back(parent_name);
      }
    }
  }

  return parents;
}

/**
 * Loop over every vertex as a starting point and then traverse the graph to the end.
 * If a vertex is encountered twice, then the graph is cyclic.
 */
bool Graph::is_cyclic()
{
  for (auto& pair : m_vertices)
  {
    // Start at every possible vertex and search all children for a loop
    std::string name = pair.first;
    std::vector<std::string> visited_vertices = { name };
    
    // As long as there are children to traverse, keep going
    std::vector<std::string> vertices_to_traverse = { name };
    while (vertices_to_traverse.size() > 0)
    {
      // grab the children of each node to traverse and insert them to the list
      std::vector<std::string> children_to_traverse;
      for (size_t i = 0; i < vertices_to_traverse.size(); ++i)
      {
        std::vector<std::string> vertex_children = m_edges[vertices_to_traverse[i]];
        children_to_traverse.insert(children_to_traverse.end(), std::make_move_iterator(vertex_children.begin()), std::make_move_iterator(vertex_children.end()));
      }

      // Check for visited vertices (i.e. cycles)
      for (size_t i = 0; i < children_to_traverse.size(); ++i)
      {
        for (size_t j = 0; j < visited_vertices.size(); ++j)
        {
          if (children_to_traverse[i] == visited_vertices[j])
          {
            return true;
          }
        }
      }

      // Prepare for the next depth search
      vertices_to_traverse = children_to_traverse;
    }
  }

  // If none of the vertices belong to a cycle, then the graph contains no cycles
  return false;
}

} // namespace cldera
