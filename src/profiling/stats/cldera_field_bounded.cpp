#include "profiling/stats/cldera_field_bounded.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

#include <ekat/util/ekat_string_utils.hpp>
namespace cldera {

FieldBounded::
FieldBounded(const ekat::Comm& comm, const ekat::ParameterList& pl)
  : FieldSinglePartStat(comm,pl)
  , m_bounds(pl.get<std::vector<Real>>("Bounds"))
  , m_mask_val(m_params.get("Mask Value",0.0))

{
  // Nothing to do here
}

void FieldBounded::
compute_impl () {
  const auto dt = m_field.data_type();
  const auto rank = m_field.layout().rank();
  EKAT_REQUIRE_MSG (rank>0 && rank<=3,
      "Error! Unsupported rank in FieldBounded.\n"
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
    EKAT_ERROR_MSG ("[FieldBounded] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
  }
}

template <typename T, int N>
void FieldBounded::
do_compute_impl () {
  // Get n-dimensional stat field
  auto stat_view = m_stat_field.nd_view_nonconst<T,N>();
  Kokkos::deep_copy(stat_view, 0);

  const int part_dim = m_field.part_dim();
  const auto non_part_layout = m_field.layout().strip_dim(part_dim);
  for (int ipart = 0; ipart < m_field.nparts(); ++ipart) {
    auto fview = m_field.part_nd_view<const T,N>(ipart);
    const auto part_layout = m_field.part_layout(ipart);

    int part_size = part_layout.extent(part_dim);
    const int offset = m_field.part_offset(ipart);

    for (int i=0; i<part_size; ++i) {
      auto f = slice(fview,part_dim,i);
      auto s = slice(stat_view,part_dim,offset+i);
      if constexpr (N==1) {
        s () = m_bounds.contains(f()) ? f() : m_mask_val;
      } else if constexpr (N==2) {
        for (int j=0; j<non_part_layout.extent(0); ++i) {
          s (j) = m_bounds.contains(f(j)) ? f(j) : m_mask_val;
        }
      } else if constexpr (N==3) {
        for (int j=0; j<non_part_layout.extent(0); ++j) {
          for (int k=0; k<non_part_layout.extent(1); ++k) {
            s (j,k) = m_bounds.contains(f(j,k)) ? f(j,k) : m_mask_val;
          }
        }
      }
    }
  }
}

} // namespace cldera
