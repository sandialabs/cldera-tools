#ifndef CLDERA_MAX_FIELD_TEST_HPP
#define CLDERA_MAX_FIELD_TEST_HPP

#include "profiling/stats/cldera_field_stat.hpp"
#include "cldera_field_test.hpp"
#include "cldera_profiling_types.hpp"

#include <memory>

namespace cldera {

/*
 *  Field test class used to test whether a field is below the max bounds specified
 */

class MaxFieldTest : public FieldTest
{
public:
  MaxFieldTest(const std::string& name,
               const std::shared_ptr<const Field>& field, 
               const Real& max,
               const ekat::Comm& comm);

  // True if test passes, false if test fails
  bool test(const TimeStamp& t) override;

private:
  const Real m_max;

  std::shared_ptr<FieldStat>  m_max_stat;
};

} // namespace cldera

#endif /* CLDERA_MAX_FIELD_TEST_HPP */
