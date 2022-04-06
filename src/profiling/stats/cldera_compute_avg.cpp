#include "profiling/cldera_field.hpp"

namespace cldera {

Real compute_sum (const Field& f);

Real compute_avg (const Field& f)
{
  return compute_sum (f) / f.layout().size();
}

} // namespace cldera
