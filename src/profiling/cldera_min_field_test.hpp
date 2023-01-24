#ifndef CLDERA_MIN_FIELD_TEST_HPP
#define CLDERA_MIN_FIELD_TEST_HPP

#include "profiling/stats/cldera_field_stat.hpp"
#include "cldera_field_test.hpp"
#include "cldera_profiling_types.hpp"

#include <memory>

namespace cldera {

/*
 *  Field test class used to test whether a field is above the min bounds specified
 */

class MinFieldTest : public FieldTest
{
public:
  MinFieldTest(const std::string& name, 
               const std::shared_ptr<const Field>& field, 
               const Real& min,
               const ekat::Comm& comm);

  // True if test passes, false if test fails
  bool test(const TimeStamp& t) override;

private:
  const Real m_min;

  std::shared_ptr<FieldStat>  m_min_stat;
};

} // namespace cldera

#endif /* CLDERA_MIN_FIELD_TEST_HPP */
