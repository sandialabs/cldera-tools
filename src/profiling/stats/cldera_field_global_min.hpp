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
  void compute_impl (const Field& f, Field& stat) const override {
    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG ("[FieldGlobalMin] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template<typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    T min = std::numeric_limits<T>::max();
    for (int p=0; p<f.nparts(); ++p) {
      const auto& pl = f.part_layout(p);
      const auto& data = f.part_data<const T>(p);
      for (int i=0; i<pl.size(); ++i) {
        min = std::min(min,data[i]);
      }
    }

    m_comm.all_reduce(&min,stat.data_nonconst<T>(),1,MPI_MIN);
  }
  ekat::Comm    m_comm;
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_MIN_HPP

