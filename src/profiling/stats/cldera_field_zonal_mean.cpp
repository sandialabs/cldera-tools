#include "profiling/stats/cldera_field_zonal_mean.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

#include <algorithm>
#include <limits>

namespace cldera {

FieldZonalMean::
FieldZonalMean(const ekat::Comm& comm, const ekat::ParameterList& pl)
 : FieldStatAlongAxis(comm,pl,"ncol")
 , m_lat_bounds(pl.get<std::vector<Real>>("Latitude Bounds"))
 , m_lev_bounds(m_params.get("Level Bounds",std::vector<int>{0,std::numeric_limits<int>::max()}))
{
  // Nothing to do here
}

void FieldZonalMean::
set_aux_fields_impl ()
{
  check_aux_fields(get_aux_fields_names());
  m_lat  = m_aux_fields.at("lat");
  m_area = m_aux_fields.at("area");

  EKAT_REQUIRE_MSG (m_field.nparts() == m_lat.nparts(),
      "Error! Incompatible number of part for aux field 'lat'.\n"
      "  - stat name   : " + name() + "\n"
      "  - field name  : " + m_field.name() + "\n"
      "  - field nparts: " + std::to_string(m_field.nparts()) + "\n"
      "  - lat nparts  : " + std::to_string(m_lat.nparts()) + "\n");
  EKAT_REQUIRE_MSG (m_field.nparts() == m_area.nparts(),
      "Error! Incompatible number of part for aux field 'area'.\n"
      "  - stat name   : " + name() + "\n"
      "  - field name  : " + m_field.name() + "\n"
      "  - field nparts: " + std::to_string(m_field.nparts()) + "\n"
      "  - area nparts  : " + std::to_string(m_area.nparts()) + "\n");

  // Sanity checks
  EKAT_REQUIRE_MSG (m_lat.layout().size() == m_area.layout().size(),
      "Error! lat and area should be the same size!\n"
      "  - stat name  : " + name() + "\n"
      "  - lat layout : " + m_lat.layout().to_string() + "\n"
      "  - area layout: " + m_area.layout().to_string() + "\n");

  // Compute zonal area (the scaling factor of the zonal integral)
  m_zonal_area = 0.0;
  Real c = 0;
  Real temp, y;
  for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
    const auto& part_layout = m_area.part_layout(ipart);
    const auto& lat_part_data = m_lat.part_data<const Real>(ipart);
    const auto& area_part_data = m_area.part_data<const Real>(ipart);
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
  EKAT_REQUIRE_MSG (m_zonal_area > 0,
      "Error! Zonal area should be positive.\n"
      " - stat name : " << name() << "\n"
      " - zonal area: " << m_zonal_area << "\n");
}

void FieldZonalMean::create_stat_field ()
{
  FieldStat::create_stat_field();
  m_temp_memory.resize(size_of(m_field.data_type())*m_stat_field.layout().size());
}

void FieldZonalMean::
compute_impl ()
{
  EKAT_REQUIRE_MSG (m_aux_fields_set, "Error! lat/area fields not initialized!\n");

  const int nparts = m_field.nparts();

  const auto dt = m_field.data_type();
  const auto rank = m_field.layout().rank();
  EKAT_REQUIRE_MSG (rank>0 && rank<=3,
      "Error! Unsupported rank in FieldZonalMean.\n"
      " - field name: " + m_field.name() + "\n"
      " - field rank: " + std::to_string(rank) + "\n");

  if (dt==IntType) {
    if (rank==1) {
      do_compute_impl<int,1>();
    } else if (rank==2) {
      do_compute_impl<int,2>();
    } else {
      do_compute_impl<int,3>();
    }
  } else if (dt==RealType) {
    if (rank==1) {
      do_compute_impl<Real,1>();
    } else if (rank==2) {
      do_compute_impl<Real,2>();
    } else {
      do_compute_impl<Real,3>();
    }
  } else {
    EKAT_ERROR_MSG ("[FieldZonalMean] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
  }
}

template <typename T, int N>
void FieldZonalMean::
do_compute_impl ()
{
  // Get n-dimensional stat field
  auto stat_view = m_stat_field.nd_view_nonconst<T,N-1>();
  Kokkos::deep_copy(stat_view, 0);

  const auto& stat_dims = m_stat_field.layout().dims();

  const int col_dim = m_field.layout().idim("ncol");

  auto temp_c = view_Nd_host<T,N-1>(reinterpret_cast<T*>(m_temp_memory.data()), m_stat_field.layout().kokkos_layout());
  Kokkos::deep_copy(temp_c, 0);

  // Find, if present, the index of level dimension (in the layout stripped of ncol)
  const bool has_lev = m_stat_field.layout().has_dim("lev");
  const int lev_dim = has_lev ? m_stat_field.layout().idim("lev") : -1;

  // Small lambda to perform the Golub-Kahan summation
  T temp, y;
  auto update_stat = [&](const T& fval, T& sval, T& cval, const Real& area) {
    y = fval*area - cval;
    temp = sval + y;
    cval = (temp-sval) - y;
    sval = temp;
  };
  for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
    const auto  field_part_data = m_field.part_data<const T>(ipart);
    const auto  area_view = m_area.part_view<const Real>(ipart);
    const auto  lat_view = m_lat.part_view<const Real>(ipart);

    const auto& fpl = m_field.part_layout(ipart);
    const int ncols = fpl.extent("ncol");

    auto fview = m_field.part_nd_view<const T,N>(ipart);

    for (int icol=0; icol<ncols; ++icol) {
      const auto lat = lat_view(icol);
      if (not m_lat_bounds.contains(lat,true,true)) {
        continue;
      }
      const auto area = area_view (icol);

      // This are compile-time if statements, so they are compiled away
      if constexpr (N==1) {
        auto f_at_col = slice(fview,col_dim,icol);
        update_stat(f_at_col(),stat_view(),temp_c(),area);
      } else if constexpr (N==2) {
        auto f_at_col = slice(fview,col_dim,icol);
        for (int idim=0; idim<stat_dims[0]; ++idim) {
          if (lev_dim==0 and not m_lev_bounds.contains(idim,true,true)) {
            continue;
          }
          update_stat(f_at_col(idim),stat_view(idim),temp_c(idim),area);
        }
      } else {
        auto f_at_col = slice(fview,col_dim,icol);

        for (int idim=0; idim<stat_dims[0]; ++idim) {
          if (lev_dim==0 and not m_lev_bounds.contains(idim,true,true)) {
            continue;
          }
          for (int jdim=0; jdim<stat_dims[1]; ++jdim) {
            if (lev_dim==1 and not m_lev_bounds.contains(jdim,true,true)) {
              continue;
            }
            update_stat(f_at_col(idim,jdim),stat_view(idim,jdim),temp_c(idim,jdim),area);
          }
        }
      }
    }
  }
  // Clock MPI ops
  track_mpi_all_reduce(m_comm,stat_view.data(),m_stat_field.layout().size(),MPI_SUM,name());

  for (int i = 0; i < stat_view.size(); ++i)
    stat_view.data()[i] /= m_zonal_area;
}

} // namespace cldera
