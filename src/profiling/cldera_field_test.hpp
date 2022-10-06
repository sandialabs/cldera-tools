#ifndef CLDERA_FIELD_TEST_HPP
#define CLDERA_FIELD_TEST_HPP

#include <ekat/mpi/ekat_comm.hpp>

#include <string>

namespace cldera {

/*
 * Abstract class for field testing
 */

class FieldTest
{
public:
  FieldTest(const std::string& name) : m_name(name) {}
  virtual ~FieldTest() = default;

  // True if test passes, false if test fails
  virtual bool test() const = 0;

  const std::string& name() const { return m_name; }

private:
  const std::string m_name;
};

} // namespace cldera

#endif /* CLDERA_FIELD_TEST_HPP */
