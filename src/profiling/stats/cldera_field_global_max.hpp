#ifndef CLDERA_FIELD_GLOBAL_MAX_HPP
#define CLDERA_FIELD_GLOBAL_MAX_HPP

#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldGlobalMax : public FieldScalarStat
{
public:
  FieldGlobalMax (const ekat::Comm& comm,
                  const ekat::ParameterList& pl)
   : FieldScalarStat (comm,pl)
  { /* Nothing to do here */ }

  std::string type () const override { return "global_max"; }

protected:
  void compute_impl (const Field& f, Field& stat) const override {
    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      // WARNING: if you add support for stuff like unsigned int, the line
      //          of do_compute_impl that uses numeric_limits has to be changed!
      EKAT_ERROR_MSG ("[FieldGlobalMax] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template<typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    T max = -std::numeric_limits<T>::max();
    for (int p=0; p<f.nparts(); ++p) {
      const auto& pl = f.part_layout(p);
      const auto& data = f.part_data<const T>(p);
      for (int i=0; i<pl.size(); ++i) {
        max = std::max(max,data[i]);
      }
    }

    // Clock MPI ops
    track_mpi_all_reduce(m_comm,&max,stat.data_nonconst<T>(),1,MPI_MAX,name());
  }

  ekat::Comm    m_comm;
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_MAX_HPP
