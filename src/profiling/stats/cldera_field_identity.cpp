#include "cldera_field_identity.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"

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
do_compute_impl ()
{
  auto stat_view = m_stat_field.nd_view_nonconst<T,N>();

  const int nparts = m_field.nparts();
  const int part_dim = m_field.part_dim();
  auto non_part_layout = m_field.layout().strip_dim(part_dim);

  for (int p=0; p<nparts; ++p) {
    const auto& part_layout = m_field.part_layout(p);
    const int part_size = part_layout.dims()[part_dim];
    const int part_offset = m_field.part_offset(p);
    auto fpart_view = m_field.part_nd_view<const T,N>(p);

    for (int i=0; i<part_size; ++i) {
      if constexpr (N==1) {
        stat_view(i+part_offset) = fpart_view(i);
      } else {
        auto stat_slice = slice(stat_view,part_dim,i);
        auto f_slice    = slice(fpart_view,part_dim,i);
        for (int j=0; j<non_part_layout.extent(0); ++j) {
          if constexpr (N==2) {
            stat_slice(j) = f_slice(j);
          } else {
            for (int k=0; k<non_part_layout.extent(0); ++k) {
              stat_slice(j,k) = f_slice(j,k);
            }
          }
        }
      }
    }
  }
}

} // namespace cldera
