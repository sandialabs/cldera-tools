#ifndef CLDERA_PATHWAY_FACTORY_HPP
#define CLDERA_PATHWAY_FACTORY_HPP

#include "cldera_graph_factory.hpp"
#include "cldera_field_test_factory.hpp"
#include "cldera_pathway.hpp"
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
 * A class for building cldera::Pathway objects using an input file
 */
class PathwayFactory {
public:
  PathwayFactory()
  : m_filename("./cldera_profiling_config.yaml"),
    m_verbose(false)
  {}

  PathwayFactory(const std::string& filename, const std::map<std::string, std::shared_ptr<const Field> > fields, const bool verbose = false)
  : m_filename(filename),
    m_fields(fields),
    m_verbose(verbose)
  {}

  // Build the pathway based on the "Pathway" sublist of the input file
  std::shared_ptr<Pathway> build_pathway(std::ostream& out = std::cout) const;

private:
  // The filename associated with the graph that is to be built
  const std::string& m_filename;

  // The verbosity of the graph factory
  const bool m_verbose;

  // The fields
  std::map<std::string, std::shared_ptr<const Field> > m_fields;
};

} // namespace cldera

#endif // CLDERA_GRAPH_FACTORY_HPP
