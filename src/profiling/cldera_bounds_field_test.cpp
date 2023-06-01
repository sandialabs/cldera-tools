#include "cldera_bounds_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_field_stat.hpp"

namespace cldera {

BoundsFieldTest::BoundsFieldTest(const std::string& name, 
                                 const std::shared_ptr<const Field>& field,
                                 const Bounds<Real>& bounds,
                                 const ekat::Comm& comm) 
 : FieldTest(name, field)
 , m_bounds(bounds)
{
  m_max_stat = StatFactory::instance().create("global_max",comm,ekat::ParameterList("global_max"));
  m_max_stat->set_field(*m_field);
  m_max_stat->create_stat_field();

  m_min_stat = StatFactory::instance().create("global_min",comm,ekat::ParameterList("global_max"));
  m_min_stat->set_field(*m_field);
  m_min_stat->create_stat_field();
}

bool BoundsFieldTest::test(const TimeStamp& t)
{
  const auto field_min = m_min_stat->compute(t);
  if (field_min.data<Real>()[0] < m_bounds.min)
  {
    if(m_save_history_on_failure)
    {
      m_history.store(t, field_min.data<Real>()[0]);
    }
    return false;
  }

  const auto field_max = m_max_stat->compute(t);
  if (field_max.data<Real>()[0] > m_bounds.max)
  {
    if(m_save_history_on_failure)
    {
      m_history.store(t, field_max.data<Real>()[0]);
    }
    return false;
  }

  return true;
}

} // namespace cldera
