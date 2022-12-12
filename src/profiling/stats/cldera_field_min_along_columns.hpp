#ifndef CLDERA_FIELD_MIN_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_MIN_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldMinAlongColumns : public FieldStatAlongAxis
{
public:
  FieldMinAlongColumns (const ekat::Comm& comm)
   : FieldStatAlongAxis("ncol"),
     m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "min_along_columns"; }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
        "Error! The type cldera::Real is not capable of representing infinity.\n");

    const auto& stat_strides = compute_stat_strides(f.layout());

    auto stat_view = stat.view_nonconst<Real>();
    Kokkos::deep_copy(stat_view, std::numeric_limits<Real>::infinity());

    const int field_rank = f.layout().rank();
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& part_data = f.part_data<const Real>(ipart);
      const auto& part_layout = f.part_layout(ipart);
      const auto& part_dims = part_layout.dims();
      for (int part_index = 0; part_index < part_layout.size(); ++part_index) {
        const int stat_index = compute_stat_index(
            ipart, part_index, field_rank, field_part_dim, part_dims, stat_strides);
        stat_view(stat_index) = std::min(stat_view(stat_index), part_data[part_index]);
      }
    }

    // Since only columns are distributed, stat_view is the same size across ranks
    m_comm.all_reduce(stat_view.data(), stat_view.size(), MPI_MIN); // Use MPI_IN_PLACE
  }

  const ekat::Comm m_comm;
};

} // namespace cldera

#endif /* CLDERA_FIELD_MIN_ALONG_COLUMNS_HPP_ */
