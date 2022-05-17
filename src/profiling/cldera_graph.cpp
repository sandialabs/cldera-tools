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


void Graph::generate_dot_graph() const
{
  std::cout << "digraph graphname {" << std::endl;

  for(auto& pair : m_edges)
  {
    for (size_t i = 0; i < pair.second.size(); ++i)
    {
      std::cout << pair.first
                << " -> "
                << pair.second[i]
                << std::endl;
    }
  }
  std::cout << "}" << std::endl;
}


} // namespace cldera
