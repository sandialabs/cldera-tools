#ifndef SRC_PROFILING_STATS_CLDERA_FIELD_BOUNDING_BOX_HPP_
#define SRC_PROFILING_STATS_CLDERA_FIELD_BOUNDING_BOX_HPP_

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <memory>

namespace cldera {

class FieldBoundingBox : public FieldStat
{
public:
  FieldBoundingBox (const ekat::Comm& comm)
    : m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "bounding_box"; }

  void initialize (const std::shared_ptr<const Field>& lat, const std::shared_ptr<const Field>& lon) {
    EKAT_REQUIRE_MSG (lat->name() == "lat" && lon->name() == "lon",
        "Error! Field names are not lat and lon!\n");
    m_lat = lat;
    m_lon = lon;
  }

  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    //// TODO: currently use full layout, may want to trim this
    return FieldLayout(field_layout.dims(), field_layout.names());
  }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    EKAT_REQUIRE_MSG (m_lat != nullptr && m_lon != nullptr,
        "Error! lat/lon fields not initialized!\n");

    const int nparts = f.nparts();
    EKAT_REQUIRE_MSG (nparts == m_lat->nparts() && nparts == m_lon->nparts(),
        "Error! Field " + f.name() + " should have the same number of parts as lat/lon!\n");

    const auto& stat_strides = compute_stat_strides(f.layout());
    auto bounding_box_field = stat.view_nonconst();
    const int field_rank = f.layout().rank();
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < nparts; ++ipart) {
      const auto& field_part_data = f.part_data(ipart);
      const auto& lat_part_data = m_lat->part_data(ipart);
      const auto& lon_part_data = m_lon->part_data(ipart);
      const auto& field_part_layout = f.part_layout(ipart);
      const auto& field_part_dims = field_part_layout.dims();
      const auto& field_part_names = field_part_layout.names();
      for (int field_part_index = 0; field_part_index < field_part_layout.size(); ++field_part_index) {
        const int stat_index = compute_stat_index(
            ipart, field_part_index, field_rank, field_part_dim, field_part_dims, stat_strides);
        const int latlon_part_index = compute_latlon_part_index(field_part_index, field_rank, field_part_dims, field_part_names);
        const Real lat_val = lat_part_data[latlon_part_index];
        const Real lon_val = lon_part_data[latlon_part_index];
        if (lat_val > m_lat_bounds.min && lat_val < m_lat_bounds.max &&
            lon_val > m_lon_bounds.min && lon_val < m_lon_bounds.max)
          bounding_box_field(stat_index) = field_part_data[field_part_index];
      }
    }
  }

  inline std::vector<int> compute_stat_strides(const FieldLayout& field_layout) const
  {
    const int field_rank = field_layout.rank();
    const auto& field_dims = field_layout.dims();
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

  inline int compute_latlon_part_index(const int field_part_index, const int field_rank,
      const std::vector<int>& field_part_dims, const std::vector<std::string>& field_part_names) const
  {
    int latlon_part_index = field_part_index;
    for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
      if (field_part_names[axis] != "ncol") {
        latlon_part_index /= field_part_dims[axis];
      }
    }
    return latlon_part_index;
  }

  const ekat::Comm m_comm;
  std::shared_ptr<const Field> m_lat, m_lon;
  //// TODO: These should be input parameters
  // const Bounds m_lat_bounds = {14.0, 16.0}; // for HSW++
  // const Bounds m_lon_bounds = {119.0, 122.0}; // for HSW++
  const Bounds m_lat_bounds = {0.0, 1.0}; // for testing
  const Bounds m_lon_bounds = {0.0, 1.0}; // for testing
};

} // namespace cldera

#endif /* SRC_PROFILING_STATS_CLDERA_FIELD_BOUNDING_BOX_HPP_ */
