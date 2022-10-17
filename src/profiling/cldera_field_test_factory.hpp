#ifndef CLDERA_FIELD_TEST_FACTORY_HPP
#define CLDERA_FIELD_TEST_FACTORY_HPP

#include "cldera_field_test.hpp"
#include "cldera_min_field_test.hpp"
#include "cldera_max_field_test.hpp"
#include "cldera_bounds_field_test.hpp"
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
namespace cldera
{

/**
 * A class for building cldera::FieldTest objects using an input file
 */
class FieldTestFactory {
public:
  FieldTestFactory(const std::string& filename, 
                   const std::map<std::string, std::shared_ptr<const Field> > fields, 
                   const ekat::Comm& comm,
                   const bool verbose = false)
  : m_filename(filename),
    m_verbose(verbose),
    m_fields(fields),
    m_comm(comm)
  {}

  // Build the graph based on the "Test" sublist of the input file
  std::map<std::string, std::shared_ptr<FieldTest> > build_field_tests(std::ostream& out = std::cout);

private:
  // The filename associated with the FieldTest that is to be built
  const std::string m_filename;

  // The verbosity of the FieldTest factory
  const bool m_verbose;

  // The corresponding fields to run tests on
  std::map<std::string, std::shared_ptr<const Field> > m_fields;

  const ekat::Comm m_comm;
};

} // namespace cldera

#endif // CLDERA_FIELD_TEST_FACTORY_HPP
