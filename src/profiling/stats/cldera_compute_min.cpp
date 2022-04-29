#include "profiling/cldera_field.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

Real compute_min (const Field& f, const ekat::Comm& comm)
{
  EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
      "Error! The type cldera::Real is not capable of representing infinity.\n");

  Real min = std::numeric_limits<Real>::infinity();
  for (int p=0; p<f.nparts(); ++p) {
    const auto& pl = f.part_layout(p);
    const auto& data = f.get_part_data(p);
    for (int i=0; i<pl.size(); ++i) {
      min = std::min(min,data[i]);
    }
  }

  Real global_min;
  comm.all_reduce(&min,&global_min,1,MPI_MIN);

  return global_min;
}

} // namespace cldera
