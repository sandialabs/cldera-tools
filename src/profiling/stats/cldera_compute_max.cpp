#include "profiling/cldera_field.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

Real compute_max (const Field& f, const ekat::Comm& comm)
{
  EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
      "Error! The type cldera::Real is not capable of representing infinity.\n");

  Real max = -std::numeric_limits<Real>::infinity();
  for (int p=0; p<f.nparts(); ++p) {
    const auto& pl = f.part_layout(p);
    const auto& data = f.get_part_data(p);
    for (int i=0; i<pl.size(); ++i) {
      max = std::max(max,data[i]);
    }
  }

  Real global_max;
  comm.all_reduce(&max,&global_max,1,MPI_MAX);

  return global_max;
}

} // namespace cldera
