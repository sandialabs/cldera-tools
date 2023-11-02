#include "profiling/stats/cldera_field_bounding_box.hpp"

#include <algorithm>
#include <limits>
#include <memory>

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
void FieldBoundingBox::do_compute_impl () {
  using pair_t = Kokkos::pair<int,int>;
  using data_nd_t = typename ekat::DataND<T,N>::type;
  using strided_view_t = Kokkos::View<data_nd_t,Kokkos::LayoutStride,ekat::HostDevice,Kokkos::MemoryUnmanaged>;

  constexpr auto ALL = Kokkos::ALL();

  auto stat_view = m_stat_field.nd_view_nonconst<T,N>();

  // Determine if field has levels
  const auto& field_names = m_field.layout().names();
  const int lev_dim = std::distance(field_names.begin(), std::find(field_names.begin(), field_names.end(), "lev"));
  const bool has_lev = lev_dim<N;
  const int part_dim = m_field.part_dim();

  for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
    auto f_part_view = m_field.part_nd_view<const T,N>(ipart);
    auto lat_part_view = m_lat.part_nd_view<const Real,1>(ipart);
    auto lon_part_view = m_lon.part_nd_view<const Real,1>(ipart);

    // Helper lambda for 2d/3d cases: deduce col/lev index, and check if inside bounds
    auto in_bounds = [&](int i, int j, int k = -1) {
      auto latlon_idx = part_dim==0 ? i : (part_dim==1 ? j : k);
      auto lev_idx = lev_dim==0 ? i : (lev_dim==1 ? j : k);
      return m_lat_bounds.contains(lat_part_view(latlon_idx)) and
             m_lon_bounds.contains(lon_part_view(latlon_idx)) and
             (not has_lev or m_lev_bounds.contains(lev_idx));
    };

    const auto& part_layout = m_field.part_layout(ipart);
    const auto& dims = part_layout.dims();
    const int part_size = part_layout.dims()[part_dim];
    const int part_offset = m_field.part_offset(ipart);
    pair_t part_range(part_offset,part_offset+part_size);

    if constexpr (N==1) {
      auto stat_subview = Kokkos::subview(stat_view,part_range);
      for (int i=0; i<dims[0]; ++i) {
        if (m_lat_bounds.contains(lat_part_view(i)) and
            m_lon_bounds.contains(lon_part_view(i))) {
          stat_subview(i) = f_part_view(i);
        }
      }
    } else if constexpr (N==2) {
      auto stat_subview = part_dim==0
                        ? Kokkos::subview(stat_view,part_range,ALL)
                        : Kokkos::subview(stat_view,ALL,part_range);
      for (size_t i=0; i<stat_subview.extent(0); ++i) {
        for (size_t j=0; j<stat_subview.extent(1); ++j) {
          if (in_bounds(i,j)) {
            stat_subview(i,j) = f_part_view(i,j);
          }
        }
      }
    } else {
      // We can't use ternary op here, since the return type for the case part_dim=2
      // has LayoutStride (recall that ?: op requires both values to have the same type)
      strided_view_t stat_subview;
      switch (part_dim) {
        case 0: stat_subview = Kokkos::subview(stat_view,part_range,ALL,ALL); break;
        case 1: stat_subview = Kokkos::subview(stat_view,ALL,part_range,ALL); break;
        case 2: stat_subview = Kokkos::subview(stat_view,ALL,ALL,part_range); break;
      }
      for (size_t i=0; i<stat_subview.extent(0); ++i) {
        for (size_t j=0; j<stat_subview.extent(1); ++j) {
          for (size_t k=0; k<stat_subview.extent(2); ++k) {
            if (in_bounds(i,j,k)) {
              stat_subview(i,j,k) = f_part_view(i,j,k);
            }
          }
        }
      }
    }
  }
}

} // namespace cldera
