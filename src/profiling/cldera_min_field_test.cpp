#include "cldera_min_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_field_stat_factory.hpp"

namespace cldera {

MinFieldTest::MinFieldTest(const std::string& name, 
                           const std::shared_ptr<const Field>& field,
                           const Real& min,
                           const ekat::Comm& comm)
: FieldTest(name, field)
, m_min(min)
{
  m_min_stat = create_stat("global_min",comm);
}

bool MinFieldTest::test(const TimeStamp& t)
{
  const auto field_min = m_min_stat->compute(*m_field);
  if (field_min.data<Real>()[0] < m_min)
  {
    if(m_save_history_on_failure)
      m_history.store(t, field_min.data<Real>()[0]);
    return false;
  }

  return true;
}

} // namespace cldera
