#include "cldera_max_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_field_stat.hpp"

namespace cldera {

MaxFieldTest::MaxFieldTest(const std::string& name, 
                           const Field& field,
                           const Real& max,
                           const ekat::Comm& comm)
: FieldTest(name, field)
, m_max(max)
{
  m_max_stat = StatFactory::instance().create("global_max",comm,ekat::ParameterList("global_max"));
  m_max_stat->set_field(m_field);
  m_max_stat->create_stat_field();
}

bool MaxFieldTest::test(const TimeStamp& t)
{
  const auto field_max = m_max_stat->compute(t);
  if (field_max.data<Real>()[0] > m_max)
  {
    if(m_save_history_on_failure)
      m_history.store(t, field_max.data<Real>()[0]);
    return false;
  }

  return true;
}

} // namespace cldera
