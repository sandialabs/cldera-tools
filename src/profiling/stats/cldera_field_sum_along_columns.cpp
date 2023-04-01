#include "profiling/stats/cldera_field_sum_along_columns.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

namespace cldera {

void FieldSumAlongColumns::
compute_impl ()
{
  const auto dt = m_field.data_type();
  const auto rank = m_field.layout().rank();
  EKAT_REQUIRE_MSG (rank>0 && rank<=3,
      "Error! Unsupported rank in FieldSumAlongColumns.\n"
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
    EKAT_ERROR_MSG ("[FieldSumAlongColumns] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
  }
}

template<typename T, int N>
void FieldSumAlongColumns::
do_compute_impl ()
{
  auto stat_view = m_stat_field.nd_view_nonconst<T,N-1>();
  Kokkos::deep_copy(stat_view, 0);

  const auto& stat_dims = m_stat_field.layout().dims();
  const int col_dim = m_field.layout().idim("ncol");

  auto temp_c = view_Nd_host<T,N-1>("sum_field", m_stat_field.layout().kokkos_layout());
  Kokkos::deep_copy(temp_c, 0);

  // Small lambda to perform the Golub-Kahan summation
  T temp, y;
  auto update_stat = [&](const T& fval, T& sval, T& cval) {
    y = fval - cval;
    temp = sval + y;
    cval = (temp-sval) - y;
    sval = temp;
  };

  for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
    const auto& fpl = m_field.part_layout(ipart);
    const int ncols = fpl.extent("ncol");

    auto fview = m_field.part_nd_view<const T,N>(ipart);

    for (int icol=0; icol<ncols; ++icol) {
      if constexpr (N==1) {
        auto f_at_col = slice(fview,col_dim,icol);
        update_stat(f_at_col(),stat_view(),temp_c());
      } else if constexpr (N==2) {
        auto f_at_col = slice(fview,col_dim,icol);
        for (int idim=0; idim<stat_dims[0]; ++idim) {
          update_stat(f_at_col(idim),stat_view(idim),temp_c(idim));
        }
      } else {
        auto f_at_col = slice(fview,col_dim,icol);
        for (int idim=0; idim<stat_dims[0]; ++idim) {
          for (int jdim=0; jdim<stat_dims[1]; ++jdim) {
            update_stat(f_at_col(idim,jdim),stat_view(idim,jdim),temp_c(idim,jdim));
          }
        }
      }
    }
  }

  // Since only columns are distributed, sum_field is the same size across ranks
  // Clock MPI ops
  track_mpi_all_reduce(m_comm,stat_view.data(),stat_view.size(),MPI_SUM,name());
}

} // namespace cldera
