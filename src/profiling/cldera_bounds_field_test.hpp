#ifndef CLDERA_BOUNDS_FIELD_TEST_HPP
#define CLDERA_BOUNDS_FIELD_TEST_HPP

#include "cldera_field_test.hpp"
#include "cldera_profiling_types.hpp"

#include <memory>

namespace cldera {

class Field;

// A pair used to store min/max bounds
struct Bounds {
  Real min;
  Real max;
};

/*
 *  Field test class used to test whether a field is within the min/max bounds specified
 */

class BoundsFieldTest : public FieldTest
{
public:
  BoundsFieldTest(const std::string& name, const std::shared_ptr<const Field> field, const Bounds& bounds);

  // True if test passes, false if test fails
  bool test(const ekat::Comm& comm) const override;

private:
  const std::shared_ptr<const Field> m_field;
  const Bounds m_bounds;
};

} // namespace cldera

#endif /* CLDERA_BOUNDS_FIELD_TEST_HPP */