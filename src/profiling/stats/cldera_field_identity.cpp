#include "cldera_field_identity.hpp"

namespace cldera
{

void FieldIdentity::
compute_impl () {
  const auto dt = m_field.data_type();
  const int rank = m_field.layout().rank();
  EKAT_REQUIRE_MSG (rank<=3,
      "[FieldIdenity] Unsupported field rank.\n"
      " - field name: " + m_field.name() + "\n"
      " - field rank: " + std::to_string(rank) + "\n");

  if (dt==IntType) {
    if (rank==1) {
      do_compute_impl<int,1>();
    } else if (rank==2) {
      do_compute_impl<int,2>();
    } else if (rank==3) {
      do_compute_impl<int,3>();
    }
  } else if (dt==RealType) {
    if (rank==1) {
      do_compute_impl<Real,1>();
    } else if (rank==2) {
      do_compute_impl<Real,2>();
    } else if (rank==3) {
      do_compute_impl<Real,3>();
    }
  } else {
    EKAT_ERROR_MSG (
        "[FieldIdentity] Unrecognized/unsupported data type\n"
        " - field name: " + m_field.name() + "\n"
        " - data type : " + e2str(dt) + "\n");
  }
}

template<typename T, int N>
void FieldIdentity::
do_compute_impl () {
  using pair_t = Kokkos::pair<int,int>;
  using data_nd_t = typename ekat::DataND<T,N>::type;
  using strided_view_t = Kokkos::View<data_nd_t,Kokkos::LayoutStride,ekat::HostDevice,Kokkos::MemoryUnmanaged>;

  constexpr auto ALL = Kokkos::ALL();

  auto stat_view = m_stat_field.nd_view_nonconst<T,N>();

  const int nparts = m_field.nparts();
  const int part_dim = m_field.part_dim();
  for (int p=0; p<nparts; ++p) {
    const auto& part_layout = m_field.part_layout(p);
    const int part_size = part_layout.dims()[part_dim];
    const int part_offset = m_field.part_offset(p);
    pair_t part_range(part_offset,part_offset+part_size);

    auto fpart_view = m_field.part_nd_view<const T,N>(p);

    if constexpr (N==1) {
      auto stat_subview = Kokkos::subview(stat_view,part_range);
      for (int i=0; i<part_layout.extent(0); ++i) {
        stat_subview(i) = fpart_view(i);
      }
    } else if constexpr (N==2) {
      auto stat_subview = part_dim==0
                        ? Kokkos::subview(stat_view,part_range,ALL)
                        : Kokkos::subview(stat_view,ALL,part_range);
      for (size_t i=0; i<stat_subview.extent(0); ++i) {
        for (size_t j=0; j<stat_subview.extent(1); ++j) {
          stat_subview(i,j) = fpart_view(i,j);
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
            stat_subview(i,j,k) = fpart_view(i,j,k);
          }
        }
      }
    }
  }
}

} // namespace cldera
