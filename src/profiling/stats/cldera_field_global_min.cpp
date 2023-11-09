#include "cldera_field_global_min.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"
#include <limits>

namespace cldera {

FieldGlobalMin::
FieldGlobalMin (const ekat::Comm& comm,
                const ekat::ParameterList& pl)
 : FieldScalarStat (comm,pl)
{
  // Nothing to do here
}

void FieldGlobalMin::
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
void FieldGlobalMin::
do_compute_impl () {
  T min = std::numeric_limits<T>::max();
  for (int p=0; p<m_field.nparts(); ++p) {
    const auto& pl = m_field.part_layout(p);
    const auto& fview = m_field.part_nd_view<const T,N>(p);
    if constexpr (N==1) {
      for (size_t i=0; i<fview.extent(0); ++i) {
        min = std::min(min,fview[i]);
      }
    } else if constexpr (N==2) {
      for (size_t i=0; i<fview.extent(0); ++i) {
        for (size_t j=0; j<fview.extent(1); ++j) {
          min = std::min(min,fview(i,j));
      }}
    } else {
      for (size_t i=0; i<fview.extent(0); ++i) {
        for (size_t j=0; j<fview.extent(1); ++j) {
          for (size_t k=0; k<fview.extent(0); ++k) {
            min = std::min(min,fview(i,j,k));
      }}}
    }
  }

  // Clock MPI ops
  track_mpi_all_reduce(m_comm,&min,m_stat_field.data_nonconst<T>(),1,MPI_MIN,name());
}

} // namespace cldera
