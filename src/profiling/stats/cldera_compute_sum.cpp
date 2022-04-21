#include "profiling/cldera_field.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

Real compute_sum (const Field& f, const ekat::Comm& comm)
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

  Real global_sum;
  comm.all_reduce(&sum,&global_sum,1,MPI_SUM);

  return global_sum;
}

} // namespace cldera
