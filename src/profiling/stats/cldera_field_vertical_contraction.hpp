#ifndef CLDERA_FIELD_VERTICAL_CONTRACTION_HPP
#define CLDERA_FIELD_VERTICAL_CONTRACTION_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera {

class FieldVerticalContraction : public FieldStat {
public:

  FieldVerticalContraction (const ekat::Comm& comm,
                         const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  {
    m_lev_idx_bounds = m_params.get<std::vector<int>>("level_bounds");
    m_average = m_params.get<bool>("average",true);
  }

  FieldVerticalContraction (const ekat::Comm& comm,
                         const ekat::ParameterList& pl,
                         const Field& f)
   : FieldVerticalContraction(comm,pl)
  {
    set_field(f);
  }

  std::string type () const { return "vertical_contraction"; }

  FieldLayout stat_layout (const FieldLayout& fl) const override {
    if (fl.has_dim_name("lev")) {
      return fl.strip_dim("lev");
    } else if (fl.has_dim_name("ilev")) {
      return fl.strip_dim("ilev");
    } else {
      EKAT_ERROR_MSG (
          "Error! Input field layout does not appear to have vertical dimension.\n"
          " - input layout: " + ekat::join(fl.names(),",") + "\n");
    }
  }
protected:

  void set_field_impl (const Field& f) {
    const auto& fl = f.layout();
    if (fl.has_dim_name("lev")) {
      m_vert_dim_pos = fl.dim_idx("lev");
    } else if (fl.has_dim_name("ilev")) {
      m_vert_dim_pos = fl.dim_idx("ilev");
    } else {
      EKAT_ERROR_MSG (
          "Error! Input field layout does not appear to have vertical dimension.\n"
          " - input layout: " + ekat::join(fl.names(),",") + "\n");
    }
  }

  void compute_impl () override {
    switch (m_stat_field.data_type()) {
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
  void do_compute_impl () {
    const auto& fl = m_field.layout();
    auto den = m_lev_idx_bounds.max - m_lev_idx_bounds.min + 1;
    switch (fl.rank()) {
      case 1:
      {
        auto sdata = m_stat_field.data_nonconst<T>();
        auto fdata = m_field.data<T>();
        sdata[0] = 0;
        for (int lev=m_lev_idx_bounds.min; lev<=m_lev_idx_bounds.max; ++lev) {
          sdata[0] += fdata[lev];
        }
        if (m_average) {
          sdata[0] /= den;
        }
        break;
      }
      case 2:
      {
        auto sview = m_stat_field.nd_view_nonconst<T,1>();
        Kokkos::deep_copy(sview,0);
        auto offset = 0;
        for (int p=0; p<m_field.nparts(); ++p) {
          auto fpl = m_field.part_layout(p);
          auto spl = stat_layout(fpl);
          auto fview = m_field.part_nd_view<T,2>(p);
          for (int i=0; i<spl.dims()[0]; ++i) {
            for (int lev=m_lev_idx_bounds.min; lev<=m_lev_idx_bounds.max; ++lev) {
              if (m_vert_dim_pos==0) {
                sview(offset+i) += fview(lev,i);
              } else {
                sview(offset+i) += fview(i,lev);
              }
            }
            if (m_average) {
              sview(offset+i) /= den;
            }
          }
          offset += spl.dims()[0];
        }
        break;
      }
      case 3:
      {
        auto sview = m_stat_field.nd_view_nonconst<T,2>();
        auto part_dim = m_field.part_dim();
        auto part_dim_name = m_field.layout().names()[part_dim];
        int offset_i = 0;
        int offset_j = 0;
        for (int p=0; p<m_field.nparts(); ++p) {
          auto fpl = m_field.part_layout(p);
          auto spl = stat_layout(fpl);
          auto fview = m_field.part_nd_view<T,3>(p);
          for (int i=0; i<spl.dims()[0]; ++i) {
            for (int j=0; j<spl.dims()[1]; ++j) {
              for (int lev=m_lev_idx_bounds.min; lev<=m_lev_idx_bounds.max; ++lev) {
                if (m_vert_dim_pos==0) {
                  sview(i+offset_i,j+offset_j) += fview(lev,i,j);
                } else if (m_vert_dim_pos==1) {
                  sview(i+offset_i,j+offset_j) += fview(i,lev,j);
                } else {
                  sview(i+offset_i,j+offset_j) += fview(i,j,lev);
                }
              }
              if (m_average) {
                sview(i+offset_i,j+offset_j) /= den;
              }
            }
          }

          // Update offset into i/j dim of stat, based on which of the two is partitioned
          if (m_stat_field.layout().dim_idx(part_dim_name)==0) {
            offset_i += spl.dims()[0];
          } else {
            offset_j += spl.dims()[1];
          }
        }
        break;
      }
      default:
        EKAT_ERROR_MSG ("Error! [FieldVerticalContraction] Unsupported field rank (" + std::to_string (fl.rank()) + ")\n");
    }
  }

  Bounds<int>     m_lev_idx_bounds;
  // Bounds<Real>    m_pressure_lev_bounds;
  bool            m_average;
  int             m_vert_dim_pos;
};

} // namespace cldera

#endif // CLDERA_FIELD_VERTICAL_CONTRACTION_HPP
