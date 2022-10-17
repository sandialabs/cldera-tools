#include "cldera_max_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_field_stat_factory.hpp"

namespace cldera {

MaxFieldTest::MaxFieldTest(const std::string& name, 
                           const std::shared_ptr<const Field>& field,
                           const Real& max,
                           const ekat::Comm& comm)
: FieldTest(name, field)
, m_max(max)
{
  m_max_stat = create_stat("global_max",comm);
}

bool MaxFieldTest::test(const TimeStamp& t)
{
  const auto field_max = m_max_stat->compute(*m_field);
  if (field_max.data()[0] > m_max)
  {
    if(m_save_history_on_failure)
      m_history.store(t, field_max.data()[0]);
    return false;
  }

  return true;
}

} // namespace cldera
