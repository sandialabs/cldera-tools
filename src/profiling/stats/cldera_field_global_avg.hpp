#ifndef CLDERA_FIELD_GLOBAL_AVG_HPP
#define CLDERA_FIELD_GLOBAL_AVG_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldGlobalAvg : public FieldGlobalSum
{
public:
  FieldGlobalAvg (const ekat::Comm& comm)
   : FieldGlobalSum (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "global_avg"; }

protected:
  // NOTE: unlike global max/min/sum, we don't support IntType,
  //       so no need for extra template function
  void compute_impl (const Field& f, Field& stat) const override {
    // Sum
    FieldGlobalSum::compute_impl(f,stat);

    // Divide by size
    long long size = f.layout().size();
    long long global_size;
    m_comm.all_reduce(&size,&global_size,1,MPI_SUM);

    stat.data_nonconst<Real>()[0] /= global_size;
  }
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_AVG_HPP

