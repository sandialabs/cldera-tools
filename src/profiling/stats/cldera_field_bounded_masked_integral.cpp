#include "cldera_field_bounded_masked_integral.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/ekat_assert.hpp>

namespace cldera {

FieldBoundedMaskedIntegral::
FieldBoundedMaskedIntegral (const ekat::Comm& comm,
                    const ekat::ParameterList& pl)
 : FieldMaskedIntegral(comm,pl)
{
  m_has_bounds = m_params.isParameter("valid_bounds");
  if (m_has_bounds) {
    auto vec = m_params.get<std::vector<Real>>("valid_bounds");
    m_bounds.min = vec[0];
    m_bounds.max = vec[1];
  }
}

void FieldBoundedMaskedIntegral::
compute_impl () {
  if (not m_has_bounds) {
    FieldMaskedIntegral::compute_impl ();
    return;
  }
  
  const auto dt   = m_field.data_type();
  const int  rank = m_field.layout().rank();
  if (dt==DataType::RealType) {
    switch (rank) {
      case 1: return do_compute_impl<Real,1>();
    }
  } else if (dt==DataType::IntType) {
    switch (rank) {
      case 1: return do_compute_impl<int,1>();
    }
  } else {
    EKAT_ERROR_MSG ("[FieldBoundedMaskedIntegral] Error! Unexpected/unsupported field data type.\n"
        " - field name: " + m_field.name() + "\n"
        " - field data type: " + e2str(m_field.data_type()) + "\n"
        " - stat name: " + name () + "\n");
  }

  // I don't want to deal with a weight integral that differs between, say,
  // different levels of a 3d field, so for now limit to rank 1 fields.
  // If m_average=false, we *could* handle rank>1 fields, but I don't want
  // to do that until we need to. Besides, I don't think we'll ever use
  // a masked integral with m_average=false.
  EKAT_ERROR_MSG ("[FieldBoundedMaskedIntegral] Error! Unsupported field rank.\n"
      " - field name: " + m_field.name() + "\n"
      " - field data type: " + e2str(m_field.data_type()) + "\n"
      " - stat name: " + name () + "\n");
}

template<typename T, int N>
void FieldBoundedMaskedIntegral::
do_compute_impl ()
{
  // Init stat to 0
  auto sview = m_stat_field.nd_view_nonconst<Real,N>();
  Kokkos::deep_copy(sview,0.0);

  view_1d_host<const Real> w_view;
  if (m_use_weight or m_average) {
    w_view = m_weight_field.view<const Real>();
  }

  view_1d_host<Real> w_int_view;
  if (m_average) {
    w_int_view = m_weight_integral.view_nonconst<Real>();
    Kokkos::deep_copy(w_int_view,0);
  }

  const int part_dim = m_field.part_dim();
  auto m = m_mask_field.view<int>();

  for (int p=0; p<m_field.nparts(); ++p) {
    auto fpl = m_field.part_layout(p);
    auto fview = m_field.part_nd_view<T,1>(p);
    const int part_size = fpl.extent(part_dim);
    const int part_offset = m_field.part_offset(p);

    for (int i=0; i<part_size; ++i) {
      int midx = m_mask_val_to_stat_entry.at(m(part_offset+i));

      if constexpr (N==1) {
        if (m_bounds.contains(fview(i))) {
          auto w = m_use_weight or m_average ? w_view(part_offset+i) : 1;
          sview(midx) += fview(i) * w;
          if (m_average) {
            w_int_view(midx) += w;
          }
        }
      } else {
        throw std::runtime_error("NOT IMPLEMENTED!\n");
      }
    }
  }

  // Global reduction of (weighted) integral
  track_mpi_all_reduce(m_comm,sview.data(),sview.size(),MPI_SUM,name());

  if (m_average) {
    const int mask_dim = m_field.layout().dim_idx(m_mask_field.layout().names()[0]);
    track_mpi_all_reduce(m_comm,w_int_view.data(),w_int_view.size(),MPI_SUM,name()+"_weight_integral");
    auto num_mask_ids = sview.extent(mask_dim);
    for (int i=0; i<num_mask_ids; ++i) {
      auto w = w_int_view(i);
      if constexpr (N==1) {
        sview(i) /= w;
      } else {
        throw std::runtime_error("NOT IMPLEMENTED!\n");
      }
    }
  }
}

} // namespace cldera
