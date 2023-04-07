#ifndef CLDERA_FIELD_GLOBAL_SUM_HPP
#define CLDERA_FIELD_GLOBAL_SUM_HPP

#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>
#include <ekat/ekat_parameter_list.hpp>

#include <limits>

namespace cldera {

class FieldGlobalSum : public FieldScalarStat
{
public:
  FieldGlobalSum (const ekat::Comm& comm,
                  const ekat::ParameterList& pl)
   : FieldScalarStat(comm,pl)
  { /* Nothing to do here */ }

  std::string type () const override { return "global_sum"; }

protected:
  void compute_impl (const Field& f, Field& stat) const override {
    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG ("[FieldGlobalSum] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }
  
  template<typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    // Note: use Kahan summation to increase accuracy
    T sum = 0;
    T c = 0;
    T temp,y;
    for (int p=0; p<f.nparts(); ++p) {
      const auto& pl = f.part_layout(p);
      const auto& data = f.part_data<const T>(p);
      for (int i=0; i<pl.size(); ++i) {
        y = data[i] - c;
        temp = sum + y;
        c = (temp - sum) - y;
        sum = temp;
      }
    }

    // Clock MPI ops
    track_mpi_all_reduce(m_comm,&sum,stat.data_nonconst<T>(),1,MPI_SUM,name());
  }
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_SUM_HPP
