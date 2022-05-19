#ifndef CLDERA_PROFILING_TEST_MANAGER_HPP
#define CLDERA_PROFILING_TEST_MANAGER_HPP

#include <ekat/mpi/ekat_comm.hpp>

#include <map>
#include <memory>
#include <string>

namespace cldera {

class FieldTest;

/*
 * A class to manage tests for field data
 */

class ProfilingTestManager
{
public:
  ProfilingTestManager() = default;

  void add_field_test(const std::shared_ptr<FieldTest> field_test);

  bool has_field_test(const std::string& field_test_name) const;

  // True if test passes, false if test fails
  bool run_field_test(const std::string& field_test_name, const ekat::Comm& comm) const;

private:
  std::map<std::string, std::shared_ptr<FieldTest>> m_field_tests;
};

} // namespace cldera

#endif /* CLDERA_PROFILING_TEST_MANAGER_HPP */
