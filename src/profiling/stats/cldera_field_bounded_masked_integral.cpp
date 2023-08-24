#include "cldera_field_bounded_masked_integral.hpp"

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
  switch (m_field.data_type()) {
    case RealType:
      do_compute_impl<Real>();
      break;
    case IntType:
      do_compute_impl<int>();
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unexpected/unsupported field data type.\n"
          " - field name: " + m_field.name() + "\n"
          " - field data type: " + e2str(m_field.data_type()) + "\n"
          " - stat name: " + name () + "\n");
  }
}

template<typename T>
void FieldBoundedMaskedIntegral::
do_compute_impl ()
{
  if (not m_has_bounds) {
    FieldMaskedIntegral::do_compute_impl<T> ();
    return;
  }
  
  const auto& fl = m_field.layout();
  // I don't want to deal with a weight integral that differs between, say,
  // different levels of a 3d field, so for now limit to rank 1 fields.
  // If m_average=false, we *could* handle rank>1 fields, but I don't want
  // to do that until we need to. Besides, I don't think we'll ever use
  // a masked integral with m_average=false.
  EKAT_REQUIRE_MSG (fl.rank()==1,
      "Error! BoundedMaskedIntegral only works for rank-1 fields.\n"
      "  - field name: " + m_field.name() + "\n");

  // Init stat to 0
  Kokkos::deep_copy(m_stat_field.view_nonconst<Real>(),0.0);

  int midx;
  view_1d_host<const Real> w_view;
  if (m_use_weight) {
    w_view = m_weight_field.view<const Real>();
  }

  view_1d_host<Real> w_int_view;
  if (m_average) {
    w_int_view = m_weight_integral.view_nonconst<Real>();
    Kokkos::deep_copy(w_int_view,0);
  }

  auto m = m_mask_field.view<int>();

  // I know, I said I limit myself to rank-1 fields, so what's with the switch?
  // I leave a switch in case we can generalize to rank>1 fields.
  switch (fl.rank()) {
    case 1:
    {
      // Offset of part-index into the unpartitioned field
      int offset = 0;

      auto sview = m_stat_field.nd_view_nonconst<Real,1>();

      // Loop over input field parts
      for (int p=0; p<m_field.nparts(); ++p) {
        auto fpl = m_field.part_layout(p);
        auto fview = m_field.part_nd_view<T,1>(p);

        // Loop over entries of this part
        for (int i=0; i<fpl.dims()[0]; ++i) {
          if (not m_bounds.contains(fview(i))) {
            continue;
          }
          midx = m_mask_val_to_stat_entry.at(m(offset+i));

          // Update field (weighted) integral (and weight integral)
          if (m_use_weight) {
            sview(midx) += fview(i)*w_view(i+offset);
          } else {
            sview(midx) += fview(i);
          }
          if (m_average) {
            w_int_view(midx) += w_view(i+offset);
          }
        }

        // Update offset
        offset += fpl.dims()[0];
      }

      // Global reduction of (weighted) integral
      m_comm.all_reduce(sview.data(),sview.size(),MPI_SUM);
      break;
    }
    default:
      EKAT_ERROR_MSG ("Error! [FieldBoundedMaskedIntegral] Unsupported field rank (" + std::to_string (fl.rank()) + ")\n");
  }

  if (m_average) {
    m_comm.all_reduce(w_int_view.data(),w_int_view.size(),MPI_SUM);
    switch (fl.rank()) {
      case 1:
      {
        auto sview = m_stat_field.nd_view_nonconst<Real,1>();
        for (int i=0; i<sview.extent_int(0); ++i) {
          auto w_int_val = w_int_view(i);
          if (w_int_view(i)==0) {
            sview(i) = std::numeric_limits<Real>::max();
          } else {
            sview(i) /= w_int_val;
          }
        }
        break;
      }
      default:
        EKAT_ERROR_MSG ("Error! [FieldBoundedMaskedIntegral] Unsupported field rank (" + std::to_string (fl.rank()) + ")\n");
    }
  }
}

} // namespace cldera
