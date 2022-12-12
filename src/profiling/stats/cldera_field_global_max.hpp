#ifndef CLDERA_FIELD_GLOBAL_MAX_HPP
#define CLDERA_FIELD_GLOBAL_MAX_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldGlobalMax : public FieldScalarStat
{
public:
  FieldGlobalMax (const ekat::Comm& comm)
   : m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "global_max"; }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
        "Error! The type cldera::Real is not capable of representing infinity.\n");

    Real max = -std::numeric_limits<Real>::infinity();
    for (int p=0; p<f.nparts(); ++p) {
      const auto& pl = f.part_layout(p);
      const auto& data = f.part_data<const Real>(p);
      for (int i=0; i<pl.size(); ++i) {
        max = std::max(max,data[i]);
      }
    }

    m_comm.all_reduce(&max,stat.data_nonconst<Real>(),1,MPI_MAX);
  }

  ekat::Comm    m_comm;
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_MAX_HPP
