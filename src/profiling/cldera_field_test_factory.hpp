#ifndef CLDERA_FIELD_TEST_FACTORY_HPP
#define CLDERA_FIELD_TEST_FACTORY_HPP

#include "cldera_field_test.hpp"
#include "cldera_profiling_archive.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>

#include <iostream>
#include <string>

namespace cldera
{

/**
 * A class for building cldera::FieldTest objects using an input file
 */
class FieldTestFactory {
public:
  FieldTestFactory(const ekat::ParameterList& params,
                   const ekat::Comm& comm,
                   const bool verbose = false)
  : m_params(params),
    m_verbose(verbose),
    m_comm(comm)
  {}

  // Build the graph based on the "Test" sublist of the input file
  std::map<std::string, std::shared_ptr<FieldTest> >
  build_field_tests(const ProfilingArchive& archive,
                    std::ostream& out = std::cout);

private:
  // The parameter list containing params to build the FieldTest
  const ekat::ParameterList m_params;

  // The verbosity of the FieldTest factory
  const bool m_verbose;

  const ekat::Comm m_comm;
};

} // namespace cldera

#endif // CLDERA_FIELD_TEST_FACTORY_HPP
