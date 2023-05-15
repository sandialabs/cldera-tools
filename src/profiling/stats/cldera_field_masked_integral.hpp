#ifndef CLDERA_FIELD_MASKED_INTEGRAL_HPP
#define CLDERA_FIELD_MASKED_INTEGRAL_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera {

class FieldMaskedIntegral : public FieldStat {
public:

  FieldMaskedIntegral (const ekat::Comm& comm,
                       const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  {
    m_mask_value   = m_params.get<int>("mask_value");
    m_mask_dim_name = m_params.get<std::string>("mask_dim_name");
    m_use_weight = m_params.isParameter("weight_field");
  }

  std::string type () const { return "vertical_contraction"; }

  std::vector<std::string> get_aux_fields_names () const {
    std::vector<std::string> aux_fnames;
    aux_fnames.push_back(m_params.get<std::string>("mask_field"));
    if (m_use_weight) {
      aux_fnames.push_back(m_params.get<std::string>("weight_field"));
    }
    return aux_fnames;
  }

  FieldLayout stat_layout (const FieldLayout& fl) const override {
    return fl.strip_dim(m_mask_dim_name);
  }

  // Since we may have weights, let's just always use Real for the result.
  DataType stat_data_type() const { return DataType::RealType; }
protected:

  void set_aux_fields_impl (const std::map<std::string,Field>& fields) {
    // For sure we have the mask field
    const auto& mask_name = m_params.get<std::string>("mask_field");
    const auto& m = fields.at(mask_name);
    const auto& fl = m.layout();
    EKAT_REQUIRE_MSG (fl.rank()==1,
        "Error! FieldMaskedIntegral requires a 1-dim mask field.\n"
        " - stat name: " + name() + "\n"
        " - mask field name: " + m.name() + "\n"
        " - mask field layout: " + ekat::join(fl.names(),",") + "\n");

    EKAT_REQUIRE_MSG (fl.names()[0]==m_mask_dim_name,
        "Error! Input mask field layout does not match stored mask dim name.\n"
        " - stat name: " + name() + "\n"
        " - input mask dim name: " + fl.names()[0] + "\n"
        " - stored mask dim name: " + m_mask_dim_name + "\n");

    // To avoid complicated if statements at runtime, ensure m_mask_field has one part.
    // If it has 2+ parts, create a single-part version via FieldIdentity stat.
    // WARNING: the mask is likely constant over the simulation, so this makes sense.
    //          If this is not true in the future, mask must be kept "as is",
    //          or turned into single-part field at compute time.
    if (m.nparts()==1) {
      m_mask_field = m;
    } else {
      auto& s = StatFactory::instance();
      FieldIdentity id(m_comm,ekat::ParameterList(m.name()));
      id.set_field(m);
      m_mask_field = id.compute(m_timestamp);
    }

    if (m_use_weight) {
      m_weight_field = fields.at(m_params.get<std::string>("weight_field","NONE"));
      // We assume the same weight for all the field components on each of the
      // integration points.
      EKAT_REQUIRE_MSG (m_weight_field.layout()==m_mask_field.layout(),
          "Error! Weight field layout incompatible with mask layout.\n"
          " - stat name: " + name() + "\n"
          " - weight field layout: " + ekat::join(m_weight_field.layout().names(),",") + "\n"
          " - mask field layout  : " + ekat::join(m_mask_field.layout().names(),",") + "\n");
      EKAT_REQUIRE_MSG (m_weight_field.data_type()==DataType::RealType,
          "Error! The weight field should have Real data type.\n"
          " - stat name: " + name() + "\n"
          " - weight data type: " + e2str(m_weight_field.data_type()) + "\n");
      EKAT_REQUIRE_MSG (m_weight_field.nparts()==m_field.nparts(),
          "Error! The weight field must have the same number of parts as the input field.\n"
          " - stat name: " + name() + "\n"
          " - weight field nparts: " + std::to_string(m_weight_field.nparts()) + "\n"
          " - input field nparts: " + std::to_string(m_field.nparts()) + "\n");
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
    const auto& fl = m_field.layout();
    auto m = m_mask_field.view<int>();

    // We'll use this to decide which of the input field's indices
    // is needed to index the stat field
    auto mask_dim_pos = m_field.layout().dim_idx(m_mask_field.layout().names()[0]);

    view_1d_host<const Real> w;
    switch (fl.rank()) {
      case 1:
      {
        // Offset of part-index into the unpartitioned field
        int offset = 0;

        // Init stat to 0
        auto& sdata = *m_stat_field.data_nonconst<Real>();
        sdata = 0;

        // Loop over input field parts
        for (int p=0; p<m_field.nparts(); ++p) {
          auto fpl = m_field.part_layout(p);
          auto fview = m_field.part_nd_view<T,1>(p);
          if (m_use_weight) {
            w = m_weight_field.part_view<const Real>(p);
          }

          // Loop over entries of this part
          for (int i=0; i<fpl.dims()[0]; ++i) {
            if (m(offset+i)!=m_mask_value)
              continue;

            // Update field (weighted) integral (and weight integral)
            if (m_use_weight) {
              sdata += fview(i)*w(i);
            } else {
              sdata += fview(i);
            }
          }

          // Update offset
          offset += fpl.dims()[0];
        }

        // Global reduction of (weighted) integral
        m_comm.all_reduce(&sdata,1,MPI_SUM);
        break;
      }
      case 2:
      {
        // Offsets of i/j part-indices into the unpartitioned field
        int offset_i = 0;
        int offset_j = 0;

        // Init stat to 0
        auto sview = m_stat_field.nd_view_nonconst<T,1>();
        Kokkos::deep_copy(sview,0);

        // Loop over input field parts
        for (int p=0; p<m_field.nparts(); ++p) {
          auto fpl = m_field.part_layout(p);
          auto spl = stat_layout(fpl);
          auto fview = m_field.part_nd_view<T,2>(p);
          if (m_use_weight) {
            w = m_weight_field.part_view<const Real>(p);
          }

          // Loop over part indices
          for (int i=0; i<fpl.dims()[0]; ++i) {
            // Check mask only if i index is the masked one
            if (mask_dim_pos==0 && m(i+offset_i)!=m_mask_value)
              continue;

            for (int j=0; j<fpl.dims()[1]; ++j) {
              // Check mask only if j index is the masked one
              if (mask_dim_pos==1 && m(j+offset_j)!=m_mask_value)
                continue;

              // Use non-masked index to access stat,
              // and masked index to access the weight
              if (m_use_weight) {
                if (mask_dim_pos==0) {
                  sview(j+offset_j) += w(i)*fview(i,j);
                } else {
                  sview(i+offset_i) += w(j)*fview(i,j);
                }
              } else {
                // Simply add
                if (mask_dim_pos==0) {
                  sview(j+offset_j) += fview(i,j);
                } else {
                  sview(i+offset_i) += fview(i,j);
                }
              }
            }
          }

          // Update the offset of the dimension that is partitioned
          offset_i += m_field.part_dim()==0 ? fpl.dims()[0] : 0;
          offset_j += m_field.part_dim()==1 ? fpl.dims()[1] : 0;
        }

        // Global reduction of (weighted) integral
        m_comm.all_reduce(sview.data(),sview.size(),MPI_SUM);

        break;
      }
      case 3:
      {
        // Offsets of i/j/k part-indices into the unpartitioned field
        int offset_i = 0;
        int offset_j = 0;
        int offset_k = 0;

        // Init stat to 0
        auto sview = m_stat_field.nd_view_nonconst<T,2>();
        Kokkos::deep_copy(sview,0);

        // Loop over input field parts
        for (int p=0; p<m_field.nparts(); ++p) {
          auto fpl = m_field.part_layout(p);
          auto spl = stat_layout(fpl);
          auto fview = m_field.part_nd_view<T,3>(p);
          if (m_use_weight) {
            w = m_weight_field.part_view<const Real>(p);
          }

          // Loop over part indices
          for (int i=0; i<fpl.dims()[0]; ++i) {
            // Check mask only if i index is the masked one
            if (mask_dim_pos==0 && m(i+offset_i)!=m_mask_value)
              continue;

            for (int j=0; j<fpl.dims()[1]; ++j) {
              // Check mask only if j index is the masked one
              if (mask_dim_pos==1 && m(j+offset_j)!=m_mask_value)
                continue;

              for (int k=0; k<fpl.dims()[2]; ++k) {
                // Check mask only if k index is the masked one
                if (mask_dim_pos==2 && m(k+offset_k)!=m_mask_value)
                  continue;

                // Use non-masked index to access stat,
                // and masked index to access the weight
                if (m_use_weight) {
                  if (mask_dim_pos==0) {
                    sview(j+offset_j,k+offset_k) += w(i)*fview(i,j,k);
                  } else if (mask_dim_pos==1) {
                    sview(i+offset_i,k+offset_k) += w(j)*fview(i,j,k);
                  } else {
                    sview(i+offset_i,j+offset_j) += w(k)*fview(i,j,k);
                  }
                } else {
                  // Simply add
                  if (mask_dim_pos==0) {
                    sview(j+offset_j,k+offset_k) += fview(i,j,k);
                  } else if (mask_dim_pos==1) {
                    sview(i+offset_i,k+offset_k) += fview(i,j,k);
                  } else {
                    sview(i+offset_i,j+offset_j) += fview(i,j,k);
                  }
                }
              }
            }
          }

          // Update the offset of the dimension that is partitioned
          offset_i += m_field.part_dim()==0 ? fpl.dims()[0] : 0;
          offset_j += m_field.part_dim()==1 ? fpl.dims()[1] : 0;
          offset_k += m_field.part_dim()==2 ? fpl.dims()[2] : 0;
        }

        // Global reduction of (weighted) integral
        m_comm.all_reduce(sview.data(),sview.size(),MPI_SUM);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Error! [FieldMaskedIntegral] Unsupported field rank (" + std::to_string (fl.rank()) + ")\n");
    }
  }


  // The mask field
  Field         m_mask_field;

  // The value of m_mask_field that denotes entries to be tallied
  int           m_mask_value;

  // Position of mask dim name in the field layout
  int           m_mask_dim_pos;
  
  // Store the name of the dimension along which we integrate
  std::string   m_mask_dim_name;

  // Optionally, we weigh the integrand by a weight field
  bool          m_use_weight;
  Field         m_weight_field;
};

} // namespace cldera

#endif // CLDERA_FIELD_MASKED_INTEGRAL_HPP
