#include "cldera_min_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_compute_stat.hpp"

namespace cldera {

MinFieldTest::MinFieldTest(const std::string& name, const std::shared_ptr<const Field> field,
    const Real& min) : FieldTest(name, field), m_min(min)
{
}

bool MinFieldTest::test(const ekat::Comm& comm, const TimeStamp& t)
{
  const Real field_min = compute_min(*m_field, comm);
  if (field_min < m_min)
  {
    if(m_save_history_on_failure)
      m_history.store(t, field_min);
    return false;
  }

  return true;
}

} // namespace cldera
