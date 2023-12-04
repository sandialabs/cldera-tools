#ifndef CLDERA_UTILS_HPP
#define CLDERA_UTILS_HPP

#include <ekat/kokkos/ekat_kokkos_meta.hpp>
#include <ekat/ekat_type_traits.hpp>
#include <ekat/std_meta/ekat_std_type_traits.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera {

// Take a hyperslice of an Nd view. This will *always*
// return a LayoutStride view
template<typename VT>
struct SVHelper
{
  static constexpr int N = VT::Rank;
  static_assert (N>0 && N<=3, "Only rank 1, 2, and 3 implemented so far.\n");

  using value_t = typename VT::traits::value_type;
  using exec_space = typename VT::traits::execution_space;
  using sv_t = ekat::Unmanaged<Kokkos::View<typename ekat::DataND<value_t,N-1>::type,
             Kokkos::LayoutStride,exec_space>>;

  static sv_t sv (const VT& v, int idim, int k) {
    EKAT_REQUIRE_MSG (idim>=0 && idim<N,
        "Error! Subview dim out of bounds.\n"
        " - view rank: " << N << "\n"
        " - input dim: " << idim << "\n");
    EKAT_REQUIRE_MSG (k>=0 && k<v.extent_int(idim),
        "Error! Subview index out of bounds.\n"
        " - view rank: " << N << "\n"
        " - dim extent: " << v.extent(idim) << "\n"
        " - dim index: " << k << "\n");
    return sv_impl<N>(v,idim,k);
  }

  // Impl functions use SFINAE b/c Kokkos::subview checks the number of input
  // args at compile time. If/when we require C++17, we can use constexpr if,
  // which does not require false branches to compile.
  template<int R>
  typename std::enable_if<R==1,sv_t>::type
  static sv_impl (const VT& v, int /* idim */,int k) {
    return Kokkos::subview(v,k);
  }
  template<int R>
  typename std::enable_if<R==2,sv_t>::type
  static sv_impl (const VT& v, int idim,int k) {
    if (idim==0) {
      return Kokkos::subview(v,k,Kokkos::ALL);
    } else if (idim==1) {
      return Kokkos::subview(v,Kokkos::ALL,k);
    } else {
      std::abort();
    }
  }
  template<int R>
  typename std::enable_if<R==3,sv_t>::type
  static sv_impl (const VT& v, int idim,int k) {
    if (idim==0) {
      return Kokkos::subview(v,k,Kokkos::ALL,Kokkos::ALL);
    } else if (idim==1) {
      return Kokkos::subview(v,Kokkos::ALL,k,Kokkos::ALL);
    } else if (idim==2) {
      return Kokkos::subview(v,Kokkos::ALL,Kokkos::ALL,k);
    } else {
      std::abort();
    }
  }
};

template<typename VT>
typename SVHelper<VT>::sv_t
slice (const VT& v, int idim, int k)
{
  return SVHelper<VT>::sv(v,idim,k);  
}

} // namespace cldera

#endif // CLDERA_UTILS_HPP
