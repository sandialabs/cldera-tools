#ifndef SRC_PROFILING_STATS_CLDERA_FIELD_BOUNDING_BOX_HPP_
#define SRC_PROFILING_STATS_CLDERA_FIELD_BOUNDING_BOX_HPP_

#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/stats/cldera_field_stat_utils.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <algorithm>
#include <limits>
#include <memory>

namespace cldera {

class FieldBoundingBox : public FieldSinglePartStat
{
public:
  FieldBoundingBox (const ekat::Comm& comm, const ekat::ParameterList& pl)
    : m_lat_bounds(Bounds{pl.get<std::vector<Real>>("Latitude Bounds").at(0), pl.get<std::vector<Real>>("Latitude Bounds").at(1)})
    , m_lon_bounds(Bounds{pl.get<std::vector<Real>>("Longitude Bounds").at(0), pl.get<std::vector<Real>>("Longitude Bounds").at(1)})
    , m_lev_bounds(pl.isParameter("Level Bounds") ?
        Bounds{pl.get<std::vector<Real>>("Level Bounds").at(0), pl.get<std::vector<Real>>("Level Bounds").at(1)} :
        Bounds{0.0, std::numeric_limits<Real>::max()})
    , m_mask_val(pl.isParameter("Mask Value") ? pl.get<Real>("Mask Value") : 0.0)
    , m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "bounding_box"; }

  void initialize (const std::shared_ptr<const Field>& lat, const std::shared_ptr<const Field>& lon) {
    EKAT_REQUIRE_MSG (lat->name() == "lat" && lon->name() == "lon",
        "Error! Field names are not lat, lon!\n");
    m_lat = lat;
    m_lon = lon;
  }

protected:
  void compute_impl (const Field& f, Field& stat) const override {
    EKAT_REQUIRE_MSG (m_lat != nullptr && m_lon != nullptr,
        "Error! lat/lon fields not initialized!\n");

    const int nparts = f.nparts();
    EKAT_REQUIRE_MSG (nparts == m_lat->nparts() && nparts == m_lon->nparts(),
        "Error! Field " + f.name() + " should have the same number of parts as lat/lon!\n");

    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG ("[FieldBoundingBox] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template <typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    const auto& stat_strides = compute_stat_strides(f.layout());
    auto bounding_box_field = stat.view_nonconst<T>();

    // Determine if field has levels
    const int field_rank = f.layout().rank();
    const auto& field_names = f.layout().names();
    const int field_level_dim = std::distance(field_names.begin(), std::find(field_names.begin(), field_names.end(), "lev"));
    const bool has_lev = not (field_level_dim == field_rank);

    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& field_part_data = f.part_data<const T>(ipart);
      const auto& lat_part_data = m_lat->part_data<const Real>(ipart);
      const auto& lon_part_data = m_lon->part_data<const Real>(ipart);
      const auto& field_part_layout = f.part_layout(ipart);
      const auto& field_part_dims = field_part_layout.dims();
      for (int field_part_index = 0; field_part_index < field_part_layout.size(); ++field_part_index) {
        const int stat_index = compute_stat_index(
            ipart, field_part_index, field_rank, field_part_dim, field_part_dims, stat_strides);
        if (has_lev) {
          const double level_index = compute_field_dim_index(field_part_index, field_rank, field_level_dim, field_part_dims);
          if (level_index < m_lev_bounds.min || level_index > m_lev_bounds.max) {
            bounding_box_field(stat_index) = m_mask_val;
            continue;
          }
        }
        const int geo_part_index = compute_field_dim_index(field_part_index, field_rank, field_part_dim, field_part_dims);
        const Real lat_val = lat_part_data[geo_part_index];
        const Real lon_val = lon_part_data[geo_part_index];
        if (lat_val < m_lat_bounds.min || lat_val > m_lat_bounds.max ||
            lon_val < m_lon_bounds.min || lon_val > m_lon_bounds.max) {
          bounding_box_field(stat_index) = m_mask_val;
          continue;
        }
        bounding_box_field(stat_index) = field_part_data[field_part_index];
      }
    }
  }

  const Bounds m_lat_bounds, m_lon_bounds, m_lev_bounds;
  const Real m_mask_val;
  const ekat::Comm m_comm;
  std::shared_ptr<const Field> m_lat, m_lon;
};

} // namespace cldera

#endif /* SRC_PROFILING_STATS_CLDERA_FIELD_BOUNDING_BOX_HPP_ */
