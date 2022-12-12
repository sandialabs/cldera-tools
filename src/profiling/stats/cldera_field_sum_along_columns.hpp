#ifndef CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldSumAlongColumns : public FieldStatAlongAxis
{
public:
  FieldSumAlongColumns (const ekat::Comm& comm)
   : FieldStatAlongAxis("ncol"),
     m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "sum_along_columns"; }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    const auto& stat_strides = compute_stat_strides(f.layout());

    auto temp_c = view_1d_host<Real>("sum_field", stat.layout().size());
    auto stat_view = stat.view_nonconst<Real>();
    Kokkos::deep_copy(temp_c, 0.0);
    Kokkos::deep_copy(stat_view, 0.0);

    Real temp, y;
    const int field_rank = f.layout().rank();
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& part_data = f.part_data<const Real>(ipart);
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
    m_comm.all_reduce(stat_view.data(), stat.layout().size(), MPI_SUM); // Use MPI_IN_PLACE
  }

  const ekat::Comm m_comm;
};

} // namespace cldera

#endif /* CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_ */
