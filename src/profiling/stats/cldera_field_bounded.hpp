#ifndef CLDERA_FIELD_BOUNDED_HPP_
#define CLDERA_FIELD_BOUNDED_HPP_

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <memory>

namespace cldera {

class FieldBounded : public FieldStat
{
public:
  FieldBounded (const ekat::Comm& comm)
    : m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "bounded"; }

  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return FieldLayout(field_layout.dims(), field_layout.names());
  }

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
      }
    }
  }

  inline std::vector<int> compute_stat_strides(const FieldLayout& field_layout) const
  {
    const int field_rank = field_layout.rank();
    const auto& field_dims = field_layout.dims();
    const auto& field_names = field_layout.names();
    std::vector<int> stat_dims;
    std::vector<int> stat_strides(field_rank, 0);
    for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
      int stride = 1;
      for (int dim : stat_dims)
        stride *= dim;
      stat_strides[axis] = stride;
      stat_dims.push_back(field_dims[axis]);
    }
    return stat_strides;

  }

  inline int compute_stat_index(const int ipart, const int part_index, const int field_rank, const int field_part_dim,
      const std::vector<int>& part_dims, const std::vector<int>& stat_strides) const
  {
    int field_ijk, work_index = part_index, stat_index = 0;
    for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
      if (axis == field_part_dim)
        field_ijk = ipart + work_index % part_dims[axis];
      else
        field_ijk = work_index % part_dims[axis];
      stat_index += field_ijk * stat_strides[axis];
      work_index /= part_dims[axis];
    }
    return stat_index;
  }

  //// TODO: This should be an input parameter
  // const Bounds m_bounds = {0.5e-4, 1.0e-4}; // for HSW++
  const Bounds m_bounds = {12.1, 20.1}; // for testing
  const ekat::Comm m_comm;
};

} // namespace cldera

#endif /* CLDERA_FIELD_BOUNDED_HPP_ */
