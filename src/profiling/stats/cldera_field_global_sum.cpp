#include "cldera_field_global_sum.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"
#include <limits>

namespace cldera {

FieldGlobalSum::
FieldGlobalSum (const ekat::Comm& comm,
                const ekat::ParameterList& pl)
 : FieldScalarStat (comm,pl)
{
  // Nothing to do here
}

void FieldGlobalSum::
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
void FieldGlobalSum::
do_compute_impl () {
  T sum = 0;
  T c = 0;
  T temp,y;

  auto update_sum = [&](const T val) {
    y = val - c;
    temp = sum + y;
    c = (temp - sum) - y;
    sum = temp;
  };
  for (int p=0; p<m_field.nparts(); ++p) {
    const auto& pl = m_field.part_layout(p);
    const auto& fview = m_field.part_nd_view<const T,N>(p);
    if constexpr (N==1) {
      for (int i=0; i<pl.dims()[0]; ++i) {
        update_sum(fview(i));
      }
    } else if constexpr (N==2) {
      for (int i=0; i<pl.dims()[0]; ++i) {
        for (int j=0; j<pl.dims()[1]; ++j) {
          update_sum(fview(i,j));
      }}
    } else {
      for (int i=0; i<pl.dims()[0]; ++i) {
        for (int j=0; j<pl.dims()[1]; ++j) {
          for (int k=0; k<pl.dims()[2]; ++k) {
            update_sum(fview(i,j,k));
      }}}
    }
  }

  // Clock MPI ops
  track_mpi_all_reduce(m_comm,&sum,m_stat_field.data_nonconst<T>(),1,MPI_SUM,name());
}

} // namespace cldera