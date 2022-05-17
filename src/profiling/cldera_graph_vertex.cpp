#include "cldera_graph_vertex.hpp"

namespace cldera
{

void GraphVertex::print()
{
  std::cout << m_name << std::endl;
  std::cout << "Activation threshold = " << m_activation_threshold << std::endl;
  std::cout << "Deactivation threshold = " << m_deactivation_threshold << std::endl;
  std::cout << "Is active = " << m_is_active << std::endl;
}

void GraphVertex::compute_is_active(const double value)
{
  if(value >= m_activation_threshold)
    m_is_active = true;

  if(value <= m_deactivation_threshold)
    m_is_active = false;
}


} // namespace cldera
