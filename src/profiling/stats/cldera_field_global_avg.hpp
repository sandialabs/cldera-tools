#ifndef CLDERA_FIELD_GLOBAL_AVG_HPP
#define CLDERA_FIELD_GLOBAL_AVG_HPP

#include "profiling/stats/cldera_field_global_sum.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldGlobalAvg : public FieldGlobalSum
{
public:
  FieldGlobalAvg (const ekat::Comm& comm,
                  const ekat::ParameterList& pl)
   : FieldGlobalSum (comm,pl)
  { /* Nothing to do here */ }

  std::string type () const override { return "global_avg"; }

protected:
  // NOTE: unlike global max/min/sum, we don't support IntType,
  //       so no need for extra template function
  void compute_impl () override {
    // Sum
    FieldGlobalSum::compute_impl();

    // Divide by size
    long long size = m_field.layout().size();
    long long global_size;

    // Clock MPI ops
    track_mpi_all_reduce(m_comm,&size,&global_size,1,MPI_SUM,name());

    m_stat_field.data_nonconst<Real>()[0] /= global_size;
  }
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_AVG_HPP

