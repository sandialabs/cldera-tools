#include "cldera_graph_vertex.hpp"

namespace cldera
{

void GraphVertex::print(std::ostream& out)
{
  out << "GraphVertex name = " << m_name << "\n";
  out << "Is active = " << m_is_active << "\n";
}

} // namespace cldera
