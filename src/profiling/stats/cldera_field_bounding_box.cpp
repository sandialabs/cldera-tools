#include "profiling/stats/cldera_field_bounding_box.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

#include <limits>

namespace cldera {

FieldBoundingBox::
FieldBoundingBox (const ekat::Comm& comm,
                  const ekat::ParameterList& pl)
 : FieldSinglePartStat (comm,pl)
 , m_lat_bounds(pl.get<std::vector<Real>>("Latitude Bounds"))
 , m_lon_bounds(pl.get<std::vector<Real>>("Longitude Bounds"))
 , m_mask_val(m_params.get<Real>("Mask Value",0.0))
 , m_lev_bounds(m_params.get<std::vector<Real>>("Level Bounds",{0,std::numeric_limits<Real>::max()}))
{
  // Nothing to do here
}

void FieldBoundingBox::
set_aux_fields_impl ()
{
  // We really needed all the aux fields
  check_aux_fields (get_aux_fields_names());
  m_lat = m_aux_fields.at("lat");
  m_lon = m_aux_fields.at("lon");

  EKAT_REQUIRE_MSG (m_field.nparts() == m_lat.nparts(),
      "Error! Incompatible number of part for aux field 'lat'.\n"
      "  - stat name   : " + name() + "\n"
      "  - field name  : " + m_field.name() + "\n"
      "  - field nparts: " + std::to_string(m_field.nparts()) + "\n"
      "  - lat nparts  : " + std::to_string(m_lat.nparts()) + "\n");
  EKAT_REQUIRE_MSG (m_field.nparts() == m_lon.nparts(),
      "Error! Incompatible number of part for aux field 'lon'.\n"
      "  - stat name   : " + name() + "\n"
      "  - field name  : " + m_field.name() + "\n"
      "  - field nparts: " + std::to_string(m_field.nparts()) + "\n"
      "  - lon nparts  : " + std::to_string(m_lon.nparts()) + "\n");
}

void FieldBoundingBox::compute_impl ()
{
  const auto dt  = m_field.data_type();
  const int rank = m_field.layout().rank();
  EKAT_REQUIRE_MSG (rank<=3,
      "[FieldBoundingBox] Unsupported field rank.\n"
      "  - field name: " + m_field.name() + "\n"
      "  - field rank: " + std::to_string(rank) + "\n");

  if (dt==IntType) {
    switch (rank) {
      case 1: return do_compute_impl<int,1>();
      case 2: return do_compute_impl<int,2>();
      case 3: return do_compute_impl<int,3>();
    }
  } else if (dt==RealType) {
    switch (rank) {
      case 1: return do_compute_impl<Real,1>();
      case 2: return do_compute_impl<Real,2>();
      case 3: return do_compute_impl<Real,3>();
    }
  } else {
    EKAT_ERROR_MSG ("[FieldBoundingBox] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
  }
}

template <typename T,int N>
void FieldBoundingBox::do_compute_impl ()
{
  auto stat_view = m_stat_field.nd_view_nonconst<T,N>();

  const int part_dim = m_field.part_dim();

  const auto& non_part_layout = m_field.layout().strip_dim(part_dim);

  // Helper lambda for checking if 2d/3d point is inside bounds
  auto in_latlon_bounds = [&](Real lat, Real lon) {
    return m_lat_bounds.contains(lat) and
           m_lon_bounds.contains(lon);
  };

  // Determine if field has levels, and the level idx **after part dim has been stripped**
  const bool has_lev = m_field.layout().has_dim("lev");
  auto lev_idx = has_lev ? non_part_layout.dim_idx("lev") : -1;
  auto in_vert_bound = [&](int j, int k = -1) {
    switch (lev_idx) {
      case 0: return m_lev_bounds.contains(j);
      case 1: return m_lev_bounds.contains(k);
      default: return true;
    }
  };

  for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
    auto f_part_view = m_field.part_nd_view<const T,N>(ipart);
    auto lat_part_view = m_lat.part_nd_view<const Real,1>(ipart);
    auto lon_part_view = m_lon.part_nd_view<const Real,1>(ipart);

    const auto& part_layout = m_field.part_layout(ipart);
    const int part_size = part_layout.dims()[part_dim];
    const int part_offset = m_field.part_offset(ipart);

    for (int i=0; i<part_size; ++i) {
      auto lat = lat_part_view(i);
      auto lon = lon_part_view(i);
      if (not in_latlon_bounds(lat,lon)) {
        continue;
      }

      auto stat_slice = slice(stat_view,part_dim,part_offset+i);
      auto f_slice    = slice(f_part_view,part_dim,i);
      if constexpr (N==1) {
        stat_slice() = f_part_view(i);
      } else if constexpr (N==2) {
        for (int j=0; j<non_part_layout.extent(0); ++j) {
          // Note: if j is the lev dim, this check does something,
          //       otherwise it always returns true
          if (in_vert_bound(j)) {
            stat_slice(j) = f_slice(j);
          }
        }
      } else {
        for (int j=0; j<non_part_layout.extent(0); ++j) {
          for (int k=0; k<non_part_layout.extent(1); ++k) {
            if (in_vert_bound(j,k)) {
              stat_slice(j,k) = f_slice(j,k);
            }
          }
        }
      }
    }
  }
}

} // namespace cldera
