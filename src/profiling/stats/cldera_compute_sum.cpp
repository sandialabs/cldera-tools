#include "profiling/cldera_field.hpp"

namespace cldera {

Real compute_sum (const Field& f)
{
  // Note: use Kahan summation to increase accuracy
  Real sum = 0;
  Real c = 0;
  Real temp,y;
  for (int p=0; p<f.nparts(); ++p) {
    const auto& pl = f.part_layout(p);
    const auto& data = f.get_part_data(p);
    for (int i=0; i<pl.size(); ++i) {
      y = data[i] - c;
      temp = sum + y;
      c = (temp - sum) - y;
      sum = temp;
    }
  }

  return sum;
}

} // namespace cldera
