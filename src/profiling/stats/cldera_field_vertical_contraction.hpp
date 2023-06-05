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
    m_use_weight = m_params.get<std::string>("weight_field","NONE")!="NONE";

    if (m_average and m_use_weight) {
      // Use another instance of this class to compute the weight field integral
      ekat::ParameterList pl("wsum");
      pl.set("average",false);
      pl.set("level_bounds",m_lev_idx_bounds.to_vector());
      m_weight_integral_stat = std::make_shared<FieldVerticalContraction>(m_comm,pl);
    }
  }

  std::string type () const override { return "vertical_contraction"; }

  std::vector<std::string> get_aux_fields_names () const override {
    std::vector<std::string> aux_fnames;
    if (m_use_weight) {
      aux_fnames.push_back(m_params.get<std::string>("weight_field"));
    }
    return aux_fnames;
  }

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

  // For simplicity, let's just assume stat is always real
  DataType stat_data_type() const override { return DataType::RealType; }

protected:

  void set_field_impl (const Field& f) override {
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

    // Create a weight field w=1. If m_use_weight=true, it should be overwritten later
    auto lev_name = fl.names()[m_vert_dim_pos];
    auto nlev = fl.dims()[m_vert_dim_pos];
    m_weight_field = Field("weight",FieldLayout({nlev},{lev_name}),DataAccess::Copy);
    m_weight_field.commit();
    m_weight_field.deep_copy(1.0);

    // Create the 0d m_weight_integral. If user sets an actual w field, we'll recompute
    // later the correct value.
    m_weight_integral = Field("weight_sum",{},DataAccess::Copy);
    m_weight_integral.commit();
    auto& wint = m_weight_integral.data_nonconst<Real>()[0];
    if (m_average) {
      wint = m_lev_idx_bounds.max - m_lev_idx_bounds.min + 1;
    } else {
      wint = 1.0;
    }
  }

  void set_aux_fields_impl (const std::map<std::string,Field>& fields) override {
    const auto& fl = m_field.layout();
    const auto lev_dim_name = fl.names()[m_vert_dim_pos];
    const auto nlevs  = fl.dims()[m_vert_dim_pos];
    if (m_use_weight) {
      const auto& wname = m_params.get<std::string>("weight_field");
      m_weight_field = fields.at(wname);

      // Sanity checks
      EKAT_REQUIRE_MSG (m_weight_field.data_type()==RealType,
          "Error! Weight field must have real data type.\n"
          " - stat name: " + name() + "\n"
          " - weight field data type: " + e2str(m_weight_field.data_type()) + "\n");


      // Note: we allow the weight field to either be (nlev) or (nlev,ncols)
      const auto& wl = m_weight_field.layout();
      EKAT_REQUIRE_MSG (wl.rank()==1 || wl.rank()==2,
          "Error! Invalid rank for the VerticalContraction weight field.\n"
          " - stat name: " + name() + "\n"
          " - weight field layout: (" + ekat::join(wl.names(),",") + ")\n");
      EKAT_REQUIRE_MSG (fl.rank()>1 or wl.rank()==1,
          "Error! 2-dim weight field only allowed for 2+ dimensional input fields.\n"
          " - stat name: " + name() + "\n"
          " - input field layout: (" + ekat::join(fl.names(),",") + ")\n"
          " - weight field layout: (" + ekat::join(wl.names(),",") + ")\n");
      EKAT_REQUIRE_MSG (wl.has_dim_name(lev_dim_name),
          "Error! Weight field layout does not have the expected vert dim name.\n"
          " - stat name: " + name() + "\n"
          " - input field layout: (" + ekat::join(fl.names(),",") + ")\n"
          " - weight field layout: (" + ekat::join(wl.names(),",") + ")\n");

      m_weight2d = wl.rank()==2;

      // For simplicity, w must be partitioned just like m_field
      EKAT_REQUIRE_MSG (not m_weight2d or m_weight_field.nparts()==m_field.nparts(),
          "Error! 2d weight field must have same number of parts as input field.\n"
        " - stat name: " + name() + "\n"
        " - weight field layout: (" + ekat::join(wl.names(),",") + ")\n"
        " - input field nparts: " + std::to_string(m_field.nparts()) + "\n"
        " - weight field nparts: " + std::to_string(m_weight_field.nparts()) + "\n");

      if (m_weight2d and fl.rank()==3) {
        // We need to figure out which one is the weight non-lev dimension within
        // the list of dimensions of the input field
        auto non_lev_dim = wl.names()[0]==lev_dim_name ? wl.names()[1] : wl.names()[0];
        m_w_first_non_lev_dim_same_as_f = 
          m_vert_dim_pos==0 ? non_lev_dim==fl.names()[1] : non_lev_dim==fl.names()[0];

        // And we need to check that the lev and non-lev dims are ordered in the same
        // way as in the input field
        auto lev_first = fl.dim_idx(lev_dim_name)<fl.dim_idx(non_lev_dim);
        EKAT_REQUIRE_MSG ( (lev_first and wl.dim_idx(lev_dim_name)==0) or 
                           (not lev_first and wl.dim_idx(lev_dim_name)==1),
          "Error! Dimensions are ordered differently in weight and input field.\n"
          " - stat name: " + name() + "\n"
          " - input field layout: (" + ekat::join(fl.names(),",") + ")\n"
          " - weight field layout: (" + ekat::join(wl.names(),",") + ")\n");
      }

      if (m_average) {
        m_constant_weight = m_params.get("constant_weight",false);

        // Setup the weight integral stat, and, if w is constant, go ahead and compute it
        m_weight_integral_stat->set_field(m_weight_field);
        m_weight_integral_stat->create_stat_field();
        m_weight_integral = m_weight_integral_stat->get_stat_field();
        if (m_constant_weight) {
          m_weight_integral_stat->compute(m_timestamp);
        }
      }
    } else if (m_average) {
      if (m_weight2d) {
        auto lev_dim = fl.names()[m_vert_dim_pos];
        m_weight_integral = Field("wint",m_weight_field.layout().strip_dim(lev_dim),DataAccess::Copy);
      } else {
        m_weight_integral = Field("wint",FieldLayout(),DataAccess::Copy);
      }
      m_weight_integral.commit();
      m_weight_integral.deep_copy(1.0);
    }
  }

  void compute_impl () override {
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
  void do_compute_impl () {

    m_stat_field.deep_copy(0.0);
    const auto& fl = m_field.layout();

    if (m_average and not m_constant_weight) {
      m_weight_integral_stat->compute(m_timestamp);
    }
    auto wint = m_weight_integral.view<Real>();

    view_1d_host<const Real> w1d;
    view_Nd_host<const Real,2> w2d;
    if (not m_weight2d) {
      w1d = m_weight_field.view<Real>();
    }

    switch (fl.rank()) {
      case 1:
      {
        auto sdata = m_stat_field.data_nonconst<Real>();
        auto fdata = m_field.data<T>();
        for (int lev=m_lev_idx_bounds.min; lev<=m_lev_idx_bounds.max; ++lev) {
          sdata[0] += fdata[lev]*w1d[lev];
        }
        // If m_average=false, wint[0]=1
        sdata[0] /= wint[0];

        break;
      }
      case 2:
      {
        auto sview = m_stat_field.nd_view_nonconst<Real,1>();
        auto offset = 0;
        for (int p=0; p<m_field.nparts(); ++p) {
          auto fpl = m_field.part_layout(p);
          auto spl = stat_layout(fpl);
          auto fview = m_field.part_nd_view<T,2>(p);
          if (m_weight2d) {
            w2d = m_weight_field.part_nd_view<Real,2>(p);
          }
          for (int i=0; i<spl.dims()[0]; ++i) {
            for (int lev=m_lev_idx_bounds.min; lev<=m_lev_idx_bounds.max; ++lev) {
              if (m_vert_dim_pos==0) {
                auto w = m_weight2d ? w2d(lev,i) : w1d(lev);
                sview(offset+i) += fview(lev,i)*w;
              } else {
                auto w = m_weight2d ? w2d(i,lev) : w1d(lev);
                sview(offset+i) += fview(i,lev)*w;
              }
            }
            // If m_average=false, wint=1
            sview(offset+i) /= m_weight2d ? wint[offset+i] : wint[0];
          }
          offset += spl.dims()[0];
        }
        break;
      }
      case 3:
      {
        auto sview = m_stat_field.nd_view_nonconst<Real,2>();
        auto part_dim = m_field.part_dim();
        auto part_dim_name = m_field.layout().names()[part_dim];
        int offset_i = 0;
        int offset_j = 0;
        int i,j;

        // We need to figure out which non-lev dim of f (i or j) correspond
        // to the non-level dim of the weight2d (if weight is 2d).
        // E.g., with w(nlev,ncol), f1(nlev,ndim,ncol) vs f2(ncol,ndim,nlev)
        int& offset_wint = m_w_first_non_lev_dim_same_as_f ? offset_i : offset_j;
        int& wint_non_lev_i = m_w_first_non_lev_dim_same_as_f ? i : j;

        for (int p=0; p<m_field.nparts(); ++p) {
          auto fpl = m_field.part_layout(p);
          auto spl = stat_layout(fpl);
          auto fview = m_field.part_nd_view<T,3>(p);
          if (m_weight2d) {
            w2d = m_weight_field.part_nd_view<Real,2>(p);
          }
          for (i=0; i<spl.dims()[0]; ++i) {
            for (j=0; j<spl.dims()[1]; ++j) {
              for (int lev=m_lev_idx_bounds.min; lev<=m_lev_idx_bounds.max; ++lev) {
                if (m_vert_dim_pos==0) {
                  auto w = m_weight2d ? 
                    (m_w_first_non_lev_dim_same_as_f ? w2d(lev,i) : w2d(lev,j)) : w1d(lev);
                  sview(i+offset_i,j+offset_j) += fview(lev,i,j)*w;
                } else if (m_vert_dim_pos==1) {
                  auto w = m_weight2d ? 
                    (m_w_first_non_lev_dim_same_as_f ? w2d(i,lev) : w2d(lev,j)) : w1d(lev);
                  sview(i+offset_i,j+offset_j) += fview(i,lev,j)*w;
                } else {
                  auto w = m_weight2d ? 
                    (m_w_first_non_lev_dim_same_as_f ? w2d(i,lev) : w2d(j,lev)) : w1d(lev);
                  sview(i+offset_i,j+offset_j) += fview(i,j,lev)*w;
                }
              }

              // If m_average=false, wint=1
              sview(i+offset_i,j+offset_j) /=
                m_weight2d ? wint[offset_wint+wint_non_lev_i] : wint[0];
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
  int             m_vert_dim_pos;

  bool            m_use_weight;
  Field           m_weight_field;

  // Helper stat, to compute integral of weight (if needed)
  std::shared_ptr<FieldVerticalContraction> m_weight_integral_stat;

  Field           m_weight_integral;        // weight integral
  bool            m_average;                // Whether we normalize by integral of w
  bool            m_constant_weight = true; // If true, pre-compute the weight integral
  bool            m_weight2d = false;       // For 2d/3d fields, whether w is 1d or 2d

  // Only used if rank(f)=3, rank(w)=2. If true, then the non-lev dim in w corresponds
  // to the first non-lev dim in m_field, otherwise it corresponds to the second
  // non-lev dim in m_field.
  int             m_w_first_non_lev_dim_same_as_f = true;;
};

} // namespace cldera

#endif // CLDERA_FIELD_VERTICAL_CONTRACTION_HPP
