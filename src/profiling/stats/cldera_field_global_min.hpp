#ifndef CLDERA_FIELD_GLOBAL_MIN_HPP
#define CLDERA_FIELD_GLOBAL_MIN_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldGlobalMin : public FieldScalarStat
{
public:
  FieldGlobalMin (const ekat::Comm& comm)
   : m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "global_min"; }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
        "Error! The type cldera::Real is not capable of representing infinity.\n");

    Real min = std::numeric_limits<Real>::infinity();
    for (int p=0; p<f.nparts(); ++p) {
      const auto& pl = f.part_layout(p);
      const auto& data = f.part_data<const Real>(p);
      for (int i=0; i<pl.size(); ++i) {
        min = std::min(min,data[i]);
      }
    }

    m_comm.all_reduce(&min,stat.data_nonconst<Real>(),1,MPI_MIN);
  }

  ekat::Comm    m_comm;
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_MIN_HPP

