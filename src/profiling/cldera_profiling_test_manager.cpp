#include "cldera_field_test.hpp"
#include "cldera_profiling_test_manager.hpp"

#include <ekat/ekat_assert.hpp>

namespace cldera {

void ProfilingTestManager::
add_field_test(const std::shared_ptr<FieldTest> field_test)
{
  EKAT_REQUIRE_MSG(field_test != nullptr,
      "[ProfilingTestManager::add_field_test]\n"
      "  Error! Invalid field test pointer.\n");

  const auto& name = field_test->name();
  EKAT_REQUIRE_MSG(m_field_tests.find(name) == m_field_tests.end(),
      "[ProfilingTestManager::add_field_test]\n"
      "  Error! Field Test '" + name + "' was already added.\n");

  m_field_tests.emplace(name, field_test);
}

bool ProfilingTestManager::
has_field_test(const std::string& field_test_name) const
{
  return m_field_tests.find(field_test_name)!=m_field_tests.end();
}

bool ProfilingTestManager::
run_field_test(const std::string& field_test_name) const
{
  EKAT_REQUIRE_MSG(has_field_test(field_test_name),
      "[ProfilingTestManager::run_field_test] Error! Field Test '" + field_test_name + "' not found.\n");
  return m_field_tests.at(field_test_name)->test(comm,{1,1});
}

std::map<std::string, bool> ProfilingTestManager::
run_all_field_tests() const
{
  std::map<std::string, bool> results;
  for (const auto& field_test : m_field_tests)
    results[field_test.first] = field_test.second->test(comm,{1,1});
  return results;
}

} // namespace cldera
