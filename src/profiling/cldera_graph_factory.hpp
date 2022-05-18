#ifndef CLDERA_GRAPH_FACTORY_HPP
#define CLDERA_GRAPH_FACTORY_HPP

#include "cldera_graph_vertex.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "cldera_graph.hpp"

namespace cldera
{

/**
 * A class for building cldera::Graph objects using an input file
 */
class GraphFactory {
public:
  GraphFactory()
  : m_filename("./cldera_profiling_config.yaml"),
    m_verbose(false)
  {}

  GraphFactory(const std::string& filename, const bool verbose = false)
  : m_filename(filename),
    m_verbose(verbose)
  {}

  // Build the graph based on the "DAG" sublist of the input file
  std::shared_ptr<Graph> build_graph(std::ostream& out = std::cout) const;

private:
  // The filename associated with the graph that is to be built
  const std::string& m_filename;

  // The verbosity of the graph factory
  const bool m_verbose;
};

} // namespace cldera

#endif // CLDERA_GRAPH_FACTORY_HPP
