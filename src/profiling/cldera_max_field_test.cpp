#include "cldera_max_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_compute_stat.hpp"

namespace cldera {

MaxFieldTest::MaxFieldTest(const std::string& name, const std::shared_ptr<const Field> field,
    const Real& max) : FieldTest(name, field), m_max(max)
{
}

bool MaxFieldTest::test(const ekat::Comm& comm, const TimeStamp& t)
{
  const Real field_max = compute_max(*m_field, comm);
  if (field_max > m_max)
  {
    if(m_save_history_on_failure)
      m_history.store(t, field_max);
    return false;
  }

  return true;
}

} // namespace cldera
