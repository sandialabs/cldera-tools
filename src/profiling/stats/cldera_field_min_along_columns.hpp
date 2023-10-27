#ifndef CLDERA_FIELD_MIN_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_MIN_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldMinAlongColumns : public FieldStatAlongAxis
{
public:
  FieldMinAlongColumns (const ekat::Comm& comm,
                        const ekat::ParameterList& pl)
   : FieldStatAlongAxis(comm,pl,"ncol")
  { /* Nothing to do here */ }

  std::string type () const override { return "min_along_columns"; }

protected:
  void compute_impl () override {
    const auto dt = m_field.data_type();
    if (dt==IntType) {
      do_compute_impl<int>();
    } else if (dt==RealType) {
      do_compute_impl<Real>();
    } else {
      // WARNING: if you add support for stuff like unsigned int, the line
      //          of do_compute_impl that uses numeric_limits has to be changed!
      EKAT_ERROR_MSG ("[FieldMinAlongColumns] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template<typename T>
  void do_compute_impl () {
    const auto& stat_strides = compute_stat_strides(m_field.layout());

    auto stat_view = m_stat_field.view_nonconst<T>();
    Kokkos::deep_copy(stat_view, std::numeric_limits<T>::max());

    const int field_rank = m_field.layout().rank();
    const int field_part_dim = m_field.part_dim();
    for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
      const auto& part_data = m_field.part_data<const T>(ipart);
      const auto& part_layout = m_field.part_layout(ipart);
      const auto& part_dims = part_layout.dims();
      for (int part_index = 0; part_index < part_layout.size(); ++part_index) {
        const int stat_index = compute_stat_index(
            ipart, part_index, field_rank, field_part_dim, part_dims, stat_strides);
        stat_view(stat_index) = std::min(stat_view(stat_index), part_data[part_index]);
      }
    }

    // Since only columns are distributed, stat_view is the same size across ranks
    // Clock MPI ops
    track_mpi_all_reduce(m_comm,stat_view.data(),stat_view.size(),MPI_MIN,name());
  }
};

} // namespace cldera

#endif /* CLDERA_FIELD_MIN_ALONG_COLUMNS_HPP_ */
