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
  FieldBounded (const ekat::Comm& comm,
                const ekat::ParameterList& pl)
   : FieldSinglePartStat(comm,pl)
   , m_bounds(pl.get<std::vector<Real>>("Bounds"))
   , m_mask_val(m_params.get("Mask Value",0.0))
  {
    /* Nothing to do here */
  }

  std::string type () const override { return "bounded"; }

  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return field_layout;
  }

protected:
  void compute_impl () override {
    const auto dt = m_field.data_type();
    if (dt==IntType) {
      do_compute_impl<int>();
    } else if (dt==RealType) {
      do_compute_impl<Real>();
    } else {
      EKAT_ERROR_MSG ("[FieldBounded] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template <typename T>
  void do_compute_impl () {
    const auto& stat_strides = compute_stat_strides(m_field.layout());

    auto bounded_field = m_stat_field.view_nonconst<T>();

    const int field_rank = m_field.layout().rank();
    const int field_part_dim = m_field.part_dim();
    for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
      const auto& part_data = m_field.part_data<const T>(ipart);
      const auto& part_layout = m_field.part_layout(ipart);
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

  const Bounds<Real> m_bounds;
  const Real m_mask_val;
};

} // namespace cldera

#endif /* CLDERA_FIELD_BOUNDED_HPP_ */
