#ifndef CLDERA_GRAPH_VERTEX_HPP
#define CLDERA_GRAPH_VERTEX_HPP

#include "cldera_profiling_types.hpp"

#include <string>
#include <iostream>

namespace cldera
{

/**
 * A vertex class designed to store information about the state of a
 * node on a graph. A vertex is generally referred to as "active" if 
 * it has crossed above the activation threshold once and "inactive"
 * if it has crossed below the deactivation threshold or has not 
 * previously become "active".
 */
class GraphVertex {
public:
  GraphVertex (const std::string& name)
  : m_name(name),
    m_activation_threshold(0.0),
    m_deactivation_threshold(0.0),
    m_is_active(false)
  {}
  
  GraphVertex (const std::string& name, const double activation_threshold, const double deactivation_threshold)
  : m_name(name),
    m_activation_threshold(activation_threshold),
    m_deactivation_threshold(deactivation_threshold),
    m_is_active(false)
  {}
  
  // Print the class and members to std::cout (may be useful for debugging)
  void print();

  // Set is_active depending on the size of value compared to the thresholds
  void compute_is_active(const double value);

  // Setters

  // Set m_activation_threshold
  void set_activation(const double threshold) { m_activation_threshold = threshold; }

  // Set m_deactivation_threshold
  void set_deactivation(const double threshold) { m_deactivation_threshold = threshold; }

  // Set m_is_active (should probably be avoided)
  void set_is_active(const bool is_active) { m_is_active = is_active; }

  // Getters

  // Return m_name
  std::string get_name() const { return m_name; }

  // Return m_is_active
  bool get_is_active() const { return m_is_active; }

  // Return m_activation_threshold
  double get_activation_threshold() const { return m_activation_threshold; }

  // Return m_deactivation_threshold
  double get_deactivation_threshold() const { return m_deactivation_threshold; }

private:
  // The name for the object. This is const since it is used to refer to the object.
  const std::string m_name; 

  // Threshold for setting this vertex as active 
  double m_activation_threshold;

  // Threshold for setting this vertex as inactive
  double m_deactivation_threshold;

  // Boolean for specifying whether this vertex is active or not
  bool m_is_active;
};

} // namespace cldera

#endif // CLDERA_GRAPH_VERTEX_HPP
