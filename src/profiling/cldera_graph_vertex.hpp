#ifndef CLDERA_GRAPH_VERTEX_HPP
#define CLDERA_GRAPH_VERTEX_HPP

#include "cldera_profiling_types.hpp"

#include <string>
#include <iostream>

namespace cldera
{

/**
 * A vertex class designed to store information about the state of a
 * node on a graph. Currently only contains the m_is_active state.
 */
class GraphVertex {
public:
  GraphVertex (const std::string& name)
  : m_name(name),
    m_is_active(false)
  {}
  
  GraphVertex (const std::string& name, const bool is_active)
  : m_name(name),
    m_is_active(is_active)
  {}
  
  // Print the class and members to std::cout (may be useful for debugging)
  void print(std::ostream& out = std::cout);

  // Setters

  // Set m_is_active (should probably be avoided)
  void set_is_active(const bool is_active) { m_is_active = is_active; }

  // Getters

  // Return m_name
  std::string get_name() const { return m_name; }

  // Return m_is_active
  bool get_is_active() const { return m_is_active; }

private:
  // The name for the object. This is const since it is used to refer to the object.
  const std::string m_name;

  // Boolean for specifying whether this vertex is active or not
  bool m_is_active;
};

} // namespace cldera

#endif // CLDERA_GRAPH_VERTEX_HPP
