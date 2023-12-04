#include "cldera_field_sum_along_columns.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

FieldSumAlongColumns::
FieldSumAlongColumns (const ekat::Comm& comm,
                      const ekat::ParameterList& pl)
 : FieldStatAlongAxis(comm,pl,"ncol")
{
  // Nothing to do here
}

void FieldSumAlongColumns::
create_stat_field () {
  FieldStat::create_stat_field();
  m_scratch.resize (size_of(m_field.data_type()) * m_stat_field.layout().size());
}

void FieldSumAlongColumns::
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
void FieldSumAlongColumns::
do_compute_impl () {
  T* raw_scratch = reinterpret_cast<T*>(m_scratch.data());
  view_Nd_host<T,N-1> c (raw_scratch,m_stat_field.layout().kokkos_layout());
  auto stat_view = m_stat_field.nd_view_nonconst<T,N-1>();
  Kokkos::deep_copy(stat_view, 0);
  Kokkos::deep_copy(c, 0);

  T y, temp;
  auto update_sum = [&] (const T val, T& c, T& sum) {
    y = val - c;
    temp = sum + y;
    c = (temp - sum) - y;
    sum = temp;
  };

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
        update_sum(fview(icol),c(),stat_view());
      }
    } else if constexpr (N==2) {
      for (int icol=0; icol<ncols; ++icol) {
        auto f_slice = slice(fview,col_dim,icol);
        for (size_t j=0; j<stat_view.extent(0); ++j) {
          update_sum(f_slice(j),c(j),stat_view(j));
        }
      }
    } else {
      for (int icol=0; icol<ncols; ++icol) {
        auto f_slice = slice(fview,col_dim,icol);
        for (size_t j=0; j<stat_view.extent(0); ++j) {
          for (size_t k=0; k<stat_view.extent(1); ++k) {
            update_sum(f_slice(j,k),c(j,k),stat_view(j,k));
        }}
      }
    }
  }

  // Clock MPI ops
  track_mpi_all_reduce(m_comm,stat_view.data(),stat_view.size(),MPI_MAX,name());
}

} // namespace cldera
