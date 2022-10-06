#ifndef CLDERA_FIELD_GLOBAL_SUM_HPP
#define CLDERA_FIELD_GLOBAL_SUM_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldGlobalSum : public FieldScalarStat
{
public:
  FieldGlobalSum (const ekat::Comm& comm)
   : m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "global_sum"; }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    // Note: use Kahan summation to increase accuracy
    Real sum = 0;
    Real c = 0;
    Real temp,y;
    for (int p=0; p<f.nparts(); ++p) {
      const auto& pl = f.part_layout(p);
      const auto& data = f.part_data(p);
      for (int i=0; i<pl.size(); ++i) {
        y = data[i] - c;
        temp = sum + y;
        c = (temp - sum) - y;
        sum = temp;
      }
    }

    m_comm.all_reduce(&sum,stat.data_nonconst(),1,MPI_SUM);
  }

  ekat::Comm    m_comm;
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_SUM_HPP
