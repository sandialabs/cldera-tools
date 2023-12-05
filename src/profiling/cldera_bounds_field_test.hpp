#ifndef CLDERA_BOUNDS_FIELD_TEST_HPP
#define CLDERA_BOUNDS_FIELD_TEST_HPP

#include "profiling/stats/cldera_field_stat.hpp"
#include "cldera_field_test.hpp"
#include "cldera_profiling_types.hpp"

#include <memory>

namespace cldera {

/*
 *  Field test class used to test whether a field is within the min/max bounds specified
 */

class BoundsFieldTest : public FieldTest
{
public:
  BoundsFieldTest(const std::string& name,
                  const Field& field,
                  const Bounds<Real>& bounds,
                  const ekat::Comm& comm);

  // True if test passes, false if test fails
  bool test(const TimeStamp& t) override;

private:
  const Bounds<Real> m_bounds;

  std::shared_ptr<FieldStat>  m_max_stat;
  std::shared_ptr<FieldStat>  m_min_stat;
};

} // namespace cldera

#endif /* CLDERA_BOUNDS_FIELD_TEST_HPP */
