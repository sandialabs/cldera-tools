#ifndef SRC_PROFILING_STATS_CLDERA_FIELD_STAT_UTILS_HPP_
#define SRC_PROFILING_STATS_CLDERA_FIELD_STAT_UTILS_HPP_

#include <vector>

inline int compute_field_dim_index(const int field_index, const int field_rank, const int field_dim,
    const std::vector<int>& field_dims)
{
  int work_index = field_index, field_dim_index;
  for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
    if (axis == field_dim)
      field_dim_index = work_index % field_dims[field_dim];
    work_index /= field_dims[axis];
  }
  return field_dim_index;
}

inline int compute_geo_part_index(const int field_part_index, const int field_rank, const int field_part_dim,
    const std::vector<int>& field_part_dims)
{
  return compute_field_dim_index(field_part_index, field_rank, field_part_dim, field_part_dims);
}

#endif /* SRC_PROFILING_STATS_CLDERA_FIELD_STAT_UTILS_HPP_ */
