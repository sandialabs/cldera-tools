#ifndef SRC_PROFILING_STATS_CLDERA_FIELD_ZONAL_MEAN_HPP_
#define SRC_PROFILING_STATS_CLDERA_FIELD_ZONAL_MEAN_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"
#include "profiling/stats/cldera_field_stat_utils.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <algorithm>
#include <limits>
#include <memory>

namespace cldera {

class FieldZonalMean : public FieldStatAlongAxis
{
public:
  FieldZonalMean (const ekat::Comm& comm, const ekat::ParameterList& pl)
    : FieldStatAlongAxis(comm,pl,"ncol")
    , m_lat_bounds(pl.get<std::vector<Real>>("Latitude Bounds"))
    , m_lev_bounds(m_params.get("Level Bounds",std::vector<Real>{0,std::numeric_limits<Real>::max()}))
  { /* Nothing to do here */ }

  std::string type () const override { return "zonal_mean"; }

  void initialize (const std::shared_ptr<const Field>& lat, const std::shared_ptr<const Field>& area) {
    if (m_inited) {
      return;
    }

    EKAT_REQUIRE_MSG (lat->name() == "lat" && area->name() == "area",
        "Error! Field names are not lat, area!\n");
    m_lat = lat;
    m_area = area;

    const int nparts = m_area->nparts();
    EKAT_REQUIRE_MSG (nparts == m_lat->nparts(),
        "Error! area should have the same number of parts as lat!\n");

    // Compute zonal area
    EKAT_REQUIRE_MSG (m_lat->layout().size() == m_area->layout().size(),
        "Error! lat and area should be the same size!\n");
    m_zonal_area = 0.0;
    Real c = 0;
    Real temp, y;
    for (int ipart = 0; ipart < nparts; ++ipart) {
      const auto& part_layout = m_area->part_layout(ipart);
      const auto& lat_part_data = m_lat->part_data<const Real>(ipart);
      const auto& area_part_data = m_area->part_data<const Real>(ipart);
      for (int part_index = 0; part_index < part_layout.size(); ++part_index) {
        const Real lat_val = lat_part_data[part_index];
        if (lat_val > m_lat_bounds.min && lat_val < m_lat_bounds.max) {
          y = area_part_data[part_index] - c;
          temp = m_zonal_area + y;
          c = (temp - m_zonal_area) - y;
          m_zonal_area = temp;
        }
      }
    }
    track_mpi_all_reduce(m_comm,&m_zonal_area,1,MPI_SUM,name()+"_initialize");
    EKAT_REQUIRE_MSG (m_zonal_area>0,
        "Error! Zonal area is zero for stat '" + name() + "'!\n");

    m_inited = true;
  }

protected:
  void compute_impl (const Field& f, Field& stat) const override {
    EKAT_REQUIRE_MSG (m_lat != nullptr && m_area != nullptr,
        "Error! lat/area fields not initialized!\n");

    const int nparts = f.nparts();
    EKAT_REQUIRE_MSG (nparts == m_lat->nparts() && nparts == m_area->nparts(),
        "Error! Field " + f.name() + " should have the same number of parts as lat/area!\n");

    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG ("[FieldZonalMean] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template <typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    const auto& stat_strides = compute_stat_strides(f.layout());

    auto temp_c = view_1d_host<T>("temp_c", stat.layout().size());
    auto stat_view = stat.view_nonconst<T>();
    Kokkos::deep_copy(temp_c, 0);
    Kokkos::deep_copy(stat_view, 0);

    // Determine if field has levels
    const int field_rank = f.layout().rank();
    const auto& field_names = f.layout().names();
    const int field_level_dim = std::distance(field_names.begin(), std::find(field_names.begin(), field_names.end(), "lev"));
    const bool has_lev = not (field_level_dim == field_rank);

    T temp, y;
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& field_part_data = f.part_data<const T>(ipart);
      const auto& lat_part_data = m_lat->part_data<const Real>(ipart);
      const auto& area_part_data = m_area->part_data<const Real>(ipart);
      const auto& field_part_layout = f.part_layout(ipart);
      const auto& field_part_dims = field_part_layout.dims();
      for (int field_part_index = 0; field_part_index < field_part_layout.size(); ++field_part_index) {
        if (has_lev) {
          const double level_index = compute_field_dim_index(field_part_index, field_rank, field_level_dim, field_part_dims);
          if (level_index < m_lev_bounds.min || level_index > m_lev_bounds.max)
            continue;
        }
        const int geo_part_index = compute_field_dim_index(field_part_index, field_rank, field_part_dim, field_part_dims);
        const Real lat_val = lat_part_data[geo_part_index];
        if (lat_val < m_lat_bounds.min || lat_val > m_lat_bounds.max)
          continue;

        const int stat_index = compute_stat_index(
            ipart, field_part_index, field_rank, field_part_dim, field_part_dims, stat_strides);
        y = field_part_data[field_part_index] * area_part_data[geo_part_index] - temp_c(stat_index);
        temp = stat_view(stat_index) + y;
        temp_c(stat_index) = (temp - stat_view(stat_index)) - y;
        stat_view(stat_index) = temp;
      }
    }
    // Clock MPI ops
    track_mpi_all_reduce(m_comm,stat_view.data(),stat.layout().size(),MPI_SUM,name());
    for (int i = 0; i < stat_view.size(); ++i)
      stat_view(i) /= m_zonal_area;
  }

  const Bounds m_lat_bounds, m_lev_bounds;
  std::shared_ptr<const Field> m_lat, m_area;
  Real m_zonal_area = 0.0;
  bool m_inited = false;
};

} // namespace cldera

#endif /* SRC_PROFILING_STATS_CLDERA_FIELD_ZONAL_MEAN_HPP_ */
