#include "cldera_field_min_along_columns.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

#include <limits>

namespace cldera {

FieldMinAlongColumns::
FieldMinAlongColumns (const ekat::Comm& comm,
                      const ekat::ParameterList& pl)
 : FieldStatAlongAxis(comm,pl,"ncol")
{
  // Nothing to do here
}

void FieldMinAlongColumns::
compute_impl () {
  const auto dt   = m_field.data_type();
  const int  rank = m_field.layout().rank();
  if (dt==DataType::RealType) {
    switch (rank) {
      case 1: return do_compute_impl<Real,1>();
      case 2: return do_compute_impl<Real,2>();
      case 3: return do_compute_impl<Real,3>();
    }
  } else if (dt==DataType::IntType) {
    switch (rank) {
      case 1: return do_compute_impl<int,1>();
      case 2: return do_compute_impl<int,2>();
      case 3: return do_compute_impl<int,3>();
    }
  } else {
    EKAT_ERROR_MSG ("Error! Unexpected/unsupported field data type.\n"
        " - field name: " + m_field.name() + "\n"
        " - field data type: " + e2str(m_field.data_type()) + "\n"
        " - stat name: " + name () + "\n");
  }
}

template<typename T, int N>
void FieldMinAlongColumns::
do_compute_impl () {
  auto stat_view = m_stat_field.nd_view_nonconst<T,N-1>();
  Kokkos::deep_copy(stat_view, std::numeric_limits<T>::max());

  // NOTE: stat_view extents are ok, since the stat is created without
  //       padding. OTOH, the input field may be a subview of some other
  //       array, so it may contain padding. Hence, you can use stat_view's
  //       extents in your loop boundaries.
  int col_dim = m_field.layout().idim(m_axis_name);
  for (int p=0; p<m_field.nparts(); ++p) {
    const auto& pl = m_field.part_layout(p);
    const auto& fview = m_field.part_nd_view<const T,N>(p);
    const int ncols = pl.extent(m_axis_name);
    if constexpr (N==1) {
      for (int icol=0; icol<ncols; ++icol) {
        stat_view() = std::min(stat_view(),fview[icol]);
      }
    } else if constexpr (N==2) {
      for (int icol=0; icol<ncols; ++icol) {
        auto f_slice = slice(fview,col_dim,icol);
        for (size_t j=0; j<stat_view.extent(0); ++j) {
          stat_view(j) = std::min(stat_view(j),f_slice(j));
        }
      }
    } else {
      for (int icol=0; icol<ncols; ++icol) {
        auto f_slice = slice(fview,col_dim,icol);
        for (size_t j=0; j<stat_view.extent(0); ++j) {
          for (size_t k=0; k<stat_view.extent(1); ++k) {
            stat_view(j,k) = std::min(stat_view(j,k),f_slice(j,k));
        }}
      }
    }
  }

  // Clock MPI ops
  track_mpi_all_reduce(m_comm,stat_view.data(),stat_view.size(),MPI_MIN,name());
}

} // namespace cldera
