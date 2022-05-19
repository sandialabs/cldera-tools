#include "cldera_bounds_field_test.hpp"
#include "cldera_field.hpp"
#include "profiling/stats/cldera_compute_stat.hpp"

namespace cldera {

BoundsFieldTest::BoundsFieldTest(const std::string& name, const std::shared_ptr<const Field> field,
    const Bounds& bounds) : FieldTest(name), m_field(field), m_bounds(bounds)
{
}

bool BoundsFieldTest::test(const ekat::Comm& comm) const
{
  const auto field_min = compute_min(*m_field, comm);
  if (field_min < m_bounds.min)
    return false;

  const auto field_max = compute_max(*m_field, comm);
  if (field_max > m_bounds.max)
    return false;

  return true;
}

} // namespace cldera
