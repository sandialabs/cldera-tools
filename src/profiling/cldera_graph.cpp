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

} // namespace cldera
