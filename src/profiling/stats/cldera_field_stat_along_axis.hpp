#ifndef CLDERA_FIELD_STAT_ALONG_AXIS_HPP_
#define CLDERA_FIELD_STAT_ALONG_AXIS_HPP_

#include "profiling/stats/cldera_field_stat.hpp"

namespace cldera {

class FieldStatAlongAxis : public FieldStat
{
public:
  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return field_layout.strip_dim(m_axis_name);
  }

protected:
  FieldStatAlongAxis (const ekat::Comm& comm,
                      const ekat::ParameterList& params,
                      const std::string& axis_name)
   : FieldStat(comm,params)
   , m_axis_name (axis_name)
  { /* Nothing to do here */ }

  inline std::vector<int> compute_stat_strides(const FieldLayout& field_layout) const
  {
    const int field_rank = field_layout.rank();
    const auto& field_dims = field_layout.dims();
    const auto& field_names = field_layout.names();
    std::vector<int> stat_dims;
    std::vector<int> stat_strides(field_rank, 0);
    for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
      if (field_names[axis] != m_axis_name) {
        int stride = 1;
        for (int dim : stat_dims)
          stride *= dim;
        stat_strides[axis] = stride;
        stat_dims.push_back(field_dims[axis]);
      }
    }
    return stat_strides;
  }

  inline int compute_stat_index(const int ipart, const int part_index, const int field_rank, const int field_part_dim,
      const std::vector<int>& part_dims, const std::vector<int>& stat_strides) const
  {
    int field_ijk, work_index = part_index, stat_index = 0;
    for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
      if (axis == field_part_dim) {
        int part_size = part_dims[field_part_dim];
        field_ijk = ipart * part_size + work_index % part_size;
      }
      else
        field_ijk = work_index % part_dims[axis];
      stat_index += field_ijk * stat_strides[axis];
      work_index /= part_dims[axis];
    }
    return stat_index;
  }

  const std::string m_axis_name;
};

} // namespace cldera

#endif /* CLDERA_FIELD_STAT_ALONG_AXIS_HPP_ */
