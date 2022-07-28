#ifndef CLDERA_MIN_FIELD_TEST_HPP
#define CLDERA_MIN_FIELD_TEST_HPP

#include "cldera_field_test.hpp"
#include "cldera_profiling_types.hpp"

#include <memory>

namespace cldera {

class Field;

/*
 *  Field test class used to test whether a field is above the min bounds specified
 */

class MinFieldTest : public FieldTest
{
public:
  MinFieldTest(const std::string& name, const std::shared_ptr<const Field> field, const Real& min);

  // True if test passes, false if test fails
  bool test(const ekat::Comm& comm, const TimeStamp& t) override;

private:
  const Real m_min;
};

} // namespace cldera

#endif /* CLDERA_MIN_FIELD_TEST_HPP */
