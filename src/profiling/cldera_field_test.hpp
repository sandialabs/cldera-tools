#ifndef CLDERA_FIELD_TEST_HPP
#define CLDERA_FIELD_TEST_HPP

#include <ekat/mpi/ekat_comm.hpp>
#include "cldera_field.hpp"

#include <string>

namespace cldera {

/*
 * Abstract class for field testing
 */

class FieldTest
{
public:
  FieldTest(const std::string& name, const Field& field)
  : m_name(name),
    m_field(field),
    m_save_history_on_failure(false)
  {
    EKAT_REQUIRE_MSG (field.committed(),
        "Error! Input field for this test is not yet committed.\n"
        "  - test  name: " + m_name + "\n"
        "  - field name: " + field.name() + "\n");
  }

  virtual ~FieldTest() = default;

  // True if test passes, false if test fails
  virtual bool test(const TimeStamp& t) = 0;

  const std::string& name() const { return m_name; }

  const Field& get_field() const { return m_field; }

  // Set whether to save failures to history or not
  void set_save_history(const bool save_history) { m_save_history_on_failure = save_history; }

  // Get the history from the test
  History& get_test_history() { return m_history; }

protected:
  const std::string m_name;
  const Field& m_field;

  // Information related to the history of the test
  bool m_save_history_on_failure;
  History m_history;
};

} // namespace cldera

#endif /* CLDERA_FIELD_TEST_HPP */
