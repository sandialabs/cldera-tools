#include "cldera_bounds_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_field_stat_factory.hpp"

namespace cldera {

BoundsFieldTest::
BoundsFieldTest(const std::string& name,
                const std::shared_ptr<const Field>& field,
                const Bounds& bounds,
                const ekat::Comm& comm)
 : FieldTest(name)
 , m_field(field)
 , m_bounds(bounds)
{
  m_max_stat = create_stat("global_max",comm);
  m_min_stat = create_stat("global_min",comm);
}

bool BoundsFieldTest::test() const
{
  const auto field_min = m_min_stat->compute(*m_field);
  if (field_min.data()[0] < m_bounds.min)
    return false;

  const auto field_max = m_max_stat->compute(*m_field);
  if (field_max.data()[0] > m_bounds.max)
    return false;

  return true;
}

} // namespace cldera
