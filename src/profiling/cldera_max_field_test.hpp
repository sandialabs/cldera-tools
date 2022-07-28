#ifndef CLDERA_MAX_FIELD_TEST_HPP
#define CLDERA_MAX_FIELD_TEST_HPP

#include "cldera_field_test.hpp"
#include "cldera_profiling_types.hpp"

#include <memory>

namespace cldera {

class Field;

/*
 *  Field test class used to test whether a field is below the max bounds specified
 */

class MaxFieldTest : public FieldTest
{
public:
  MaxFieldTest(const std::string& name, const std::shared_ptr<const Field> field, const Real& max);

  // True if test passes, false if test fails
  bool test(const ekat::Comm& comm, const TimeStamp& t) override;

private:
  const Real m_max;
};

} // namespace cldera

#endif /* CLDERA_MAX_FIELD_TEST_HPP */
