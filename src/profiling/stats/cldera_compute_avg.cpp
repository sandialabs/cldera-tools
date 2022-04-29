#include "profiling/cldera_field.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

Real compute_sum (const Field& f,const ekat::Comm& comm);

Real compute_avg (const Field& f,const ekat::Comm& comm)
{
  auto sum = compute_sum (f,comm);
  long long size = f.layout().size();
  long long global_size;
  comm.all_reduce(&size,&global_size,1,MPI_SUM);

  return sum / global_size;
}

} // namespace cldera
