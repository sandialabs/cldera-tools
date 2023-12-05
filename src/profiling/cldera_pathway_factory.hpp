#ifndef CLDERA_PATHWAY_FACTORY_HPP
#define CLDERA_PATHWAY_FACTORY_HPP

#include "cldera_graph_factory.hpp"
#include "cldera_profiling_archive.hpp"
#include "cldera_field_test_factory.hpp"
#include "cldera_pathway.hpp"
#include "cldera_graph.hpp"

#include <ekat/ekat_parameter_list.hpp>

#include <iostream>
#include <vector>
#include <string>

namespace cldera
{

/**
 * A class for building cldera::Pathway objects using an input file
 */
class PathwayFactory {
public:
  PathwayFactory () = default;

  PathwayFactory (const ekat::ParameterList& params,
                  const ekat::Comm& comm, const bool verbose = false)
  : m_params(params),
    m_verbose(verbose),
    m_comm(comm)
  {}

  // Build the pathway based on the "Pathway" sublist of the input file
  std::shared_ptr<Pathway> build_pathway(const ProfilingArchive& archive,
                                         std::ostream& out = std::cout) const;

private:
  // The parameter list containing params to build the graph
  const ekat::ParameterList m_params;

  // The verbosity of the graph factory
  const bool m_verbose = false;

  const ekat::Comm m_comm;
};

} // namespace cldera

#endif // CLDERA_GRAPH_FACTORY_HPP
