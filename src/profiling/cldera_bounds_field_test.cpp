#include "cldera_bounds_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_field_stat_factory.hpp"

namespace cldera {

BoundsFieldTest::BoundsFieldTest(const std::string& name, 
                                 const std::shared_ptr<const Field>& field,
                                 const Bounds& bounds,
                                 const ekat::Comm& comm) 
 : FieldTest(name, field)
 , m_bounds(bounds)
{
  m_max_stat = create_stat(ekat::ParameterList("global_max"),comm);
  m_min_stat = create_stat(ekat::ParameterList("global_min"),comm);
}

bool BoundsFieldTest::test(const TimeStamp& t)
{
  const auto field_min = m_min_stat->compute(*m_field);
  if (field_min.data<Real>()[0] < m_bounds.min)
  {
    if(m_save_history_on_failure)
    {
      m_history.store(t, field_min.data<Real>()[0]);
    }
    return false;
  }

  const auto field_max = m_max_stat->compute(*m_field);
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
