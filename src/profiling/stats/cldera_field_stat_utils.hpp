#ifndef SRC_PROFILING_STATS_CLDERA_FIELD_STAT_UTILS_HPP_
#define SRC_PROFILING_STATS_CLDERA_FIELD_STAT_UTILS_HPP_

#include <vector>

inline int compute_geo_part_index(const int field_part_index, const int field_rank, const int field_part_dim,
    const std::vector<int>& field_part_dims)
{
  int work_index = field_part_index, geo_part_index;
  for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
    if (axis == field_part_dim) {
      int part_size = field_part_dims[field_part_dim];
      geo_part_index = work_index % part_size;
    }
    work_index /= field_part_dims[axis];
  }
  return geo_part_index;
}

#endif /* SRC_PROFILING_STATS_CLDERA_FIELD_STAT_UTILS_HPP_ */
