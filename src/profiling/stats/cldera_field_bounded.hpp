#ifndef CLDERA_FIELD_BOUNDED_HPP_
#define CLDERA_FIELD_BOUNDED_HPP_

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <memory>

namespace cldera {

class FieldBounded : public FieldSinglePartStat
{
public:
  FieldBounded (const ekat::Comm& comm, const ekat::ParameterList& pl)
    : m_bounds(Bounds{pl.get<std::vector<Real>>("Bounds").at(0), pl.get<std::vector<Real>>("Bounds").at(1)})
    , m_mask_val(pl.isParameter("Mask Value") ? pl.get<Real>("Mask Value") : 0.0)
    , m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "bounded"; }

protected:
  void compute_impl (const Field& f, Field& stat) const override {
    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG ("[FieldBounded] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template <typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    const auto& stat_strides = compute_stat_strides(f.layout());

    auto bounded_field = stat.view_nonconst<T>();

    const int field_rank = f.layout().rank();
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& part_data = f.part_data<const T>(ipart);
      const auto& part_layout = f.part_layout(ipart);
      const auto& part_dims = part_layout.dims();
      for (int part_index = 0; part_index < part_layout.size(); ++part_index) {
        const int stat_index = compute_stat_index(
            ipart, part_index, field_rank, field_part_dim, part_dims, stat_strides);
        const T val = part_data[part_index];
        if (val > m_bounds.min && val < m_bounds.max)
          bounded_field(stat_index) = val;
        else
          bounded_field(stat_index) = m_mask_val;
      }
    }
  }

  const Bounds m_bounds;
  const Real m_mask_val;
  const ekat::Comm m_comm;
};

} // namespace cldera

#endif /* CLDERA_FIELD_BOUNDED_HPP_ */
