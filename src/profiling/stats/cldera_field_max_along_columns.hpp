#ifndef CLDERA_FIELD_MAX_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_MAX_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldMaxAlongColumns : public FieldStatAlongAxis
{
public:
  FieldMaxAlongColumns (const ekat::Comm& comm)
   : FieldStatAlongAxis("ncol"),
     m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "max_along_columns"; }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    EKAT_REQUIRE_MSG (std::numeric_limits<Real>::has_infinity,
        "Error! The type cldera::Real is not capable of representing infinity.\n");

    const auto& stat_strides = compute_stat_strides(f.layout());

    auto max_field = view_1d_host<Real>("max_field", stat.view().size());
    Kokkos::deep_copy(max_field, -std::numeric_limits<Real>::infinity());

    const int field_rank = f.layout().rank();
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& part_data = f.part_data(ipart);
      const auto& part_layout = f.part_layout(ipart);
      const auto& part_dims = part_layout.dims();
      for (int part_index = 0; part_index < part_layout.size(); ++part_index) {
        const int stat_index = compute_stat_index(
            ipart, part_index, field_rank, field_part_dim, part_dims, stat_strides);
        max_field(stat_index) = std::max(max_field(stat_index), part_data[part_index]);
      }
    }

    // Since only columns are distributed, max_field is the same size across ranks
    m_comm.all_reduce(max_field.data(), stat.data_nonconst(), max_field.size(), MPI_MAX);
  }

  const ekat::Comm m_comm;
};

} // namespace cldera

#endif /* CLDERA_FIELD_MAX_ALONG_COLUMNS_HPP_ */
