#ifndef CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class FieldSumAlongColumns : public FieldStatAlongAxis
{
public:
  FieldSumAlongColumns (const ekat::Comm& comm,
                        const ekat::ParameterList& params)
   : FieldStatAlongAxis(comm,params,"ncol")
  { /* Nothing to do here */ }

  std::string type () const override { return "sum_along_columns"; }

protected:
  void compute_impl (const Field& f, Field& stat) const override {
    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG ("[FieldSumAlongColumns] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template<typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    const auto& stat_strides = compute_stat_strides(f.layout());

    view_1d_host<T> temp_c ("sum_field", stat.layout().size());
    auto stat_view = stat.view_nonconst<T>();
    Kokkos::deep_copy(temp_c, 0);
    Kokkos::deep_copy(stat_view, 0);

    T temp, y;
    const int field_rank = f.layout().rank();
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& part_data = f.part_data<const T>(ipart);
      const auto& part_layout = f.part_layout(ipart);
      const auto& part_dims = part_layout.dims();
      for (int part_index = 0; part_index < part_layout.size(); ++part_index) {
        const int stat_index = compute_stat_index(
            ipart, part_index, field_rank, field_part_dim, part_dims, stat_strides);
        y = part_data[part_index] - temp_c(stat_index);
        temp = stat_view(stat_index) + y;
        temp_c(stat_index) = (temp - stat_view(stat_index)) - y;
        stat_view(stat_index) = temp;
      }
    }

    // Since only columns are distributed, sum_field is the same size across ranks
    // Clock MPI ops
    track_mpi_all_reduce(m_comm,stat_view.data(),stat_view.size(),MPI_SUM,name());
  }
};

} // namespace cldera

#endif /* CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_ */
