#ifndef CLDERA_FIELD_GLOBAL_MIN_HPP
#define CLDERA_FIELD_GLOBAL_MIN_HPP

#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldGlobalMin : public FieldScalarStat
{
public:
  FieldGlobalMin (const ekat::Comm& comm,
                  const ekat::ParameterList& pl)
   : FieldScalarStat(comm,pl)
  { /* Nothing to do here */ }

  std::string type () const override { return "global_min"; }

protected:
  void compute_impl () override {
    const auto dt = m_field.data_type();
    if (dt==IntType) {
      do_compute_impl<int>();
    } else if (dt==RealType) {
      do_compute_impl<Real>();
    } else {
      EKAT_ERROR_MSG ("[FieldGlobalMin] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template<typename T>
  void do_compute_impl () {
    T min = std::numeric_limits<T>::max();
    for (int p=0; p<m_field.nparts(); ++p) {
      const auto& pl = m_field.part_layout(p);
      const auto& data = m_field.part_data<const T>(p);
      for (int i=0; i<pl.size(); ++i) {
        min = std::min(min,data[i]);
      }
    }

    // Clock MPI ops
    track_mpi_all_reduce(m_comm,&min,m_stat_field.data_nonconst<T>(),1,MPI_MIN,name());
  }
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_MIN_HPP

