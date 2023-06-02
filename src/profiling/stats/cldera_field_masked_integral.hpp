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
    m_use_weight = m_params.isParameter("weight_field");
    m_average = m_params.get<bool>("average",true);
  }

  std::string type () const { return "masked_integral"; }

  std::vector<std::string> get_aux_fields_names () const {
    std::vector<std::string> aux_fnames;
    aux_fnames.push_back(m_params.get<std::string>("mask_field"));
    if (m_use_weight) {
      aux_fnames.push_back(m_params.get<std::string>("weight_field"));
    }
    return aux_fnames;
  }

  FieldLayout stat_layout (const FieldLayout& fl) const override {
    const auto& masked_dim = m_mask_field.layout().names()[0];
    auto names = fl.names();
    auto dims  = fl.dims();
    auto pos = fl.dim_idx(masked_dim);
    dims[pos] = m_mask_val_to_stat_entry.size();
    names[pos] = "dim" + std::to_string(dims[pos]);
    return FieldLayout(dims,names);
  }

  // Since we may have weights, let's just always use Real for the result.
  DataType stat_data_type() const { return DataType::RealType; }
protected:

  void set_aux_fields_impl (const std::map<std::string,Field>& fields) {
    // For sure we have the mask field
    const auto& mask_name = m_params.get<std::string>("mask_field");
    const auto& m = fields.at(mask_name);
    const auto& fl = m.layout();
    const auto& mask_dim = fl.names()[0];
    EKAT_REQUIRE_MSG (fl.rank()==1,
        "Error! FieldMaskedIntegral requires a 1-dim mask field.\n"
        " - stat name: " + name() + "\n"
        " - mask field name: " + m.name() + "\n"
        " - mask field layout: " + ekat::join(fl.names(),",") + "\n");

    EKAT_REQUIRE_MSG (m_field.layout().has_dim_name(mask_dim),
        "Error! Input mask field layout incompatible with stored field layout.\n"
        " - stat name: " + name() + "\n"
        " - input mask layout: " + mask_dim + "\n"
        " - stored field: " + ekat::join(m_field.layout().names(),",") + "\n");

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
      id.create_stat_field();
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

    if (m_average) {
      ekat::ParameterList pl("");
      pl.set("mask_field",m_mask_field.name());
      pl.set("average",false);
      FieldMaskedIntegral w_int_stat(m_comm,pl);
      std::map<std::string,Field> aux_fields;
      aux_fields["mask"] = m_mask_field;
      if (m_use_weight) {
        w_int_stat.set_field(m_weight_field);
      } else {
        Field w("",m_mask_field.layout(),DataAccess::Copy);
        w.commit();
        Kokkos::deep_copy(w.view_nonconst<Real>(),1);
        w_int_stat.set_field(w);
      }
      w_int_stat.set_aux_fields (aux_fields);
      w_int_stat.create_stat_field();
      m_weight_integral = w_int_stat.compute(m_timestamp);
    }

    // First, gather all the mask values we have
    auto data = m_mask_field.data<int>();
    auto size = m_mask_field.layout().size();
    std::set<int> vals;
    for (int i=0; i<size; ++i) {
      vals.insert(data[i]);
    }

    // Then, map each mask value to an index in 0,...,num_mask_values-1
    for (auto v : vals) {
      m_mask_val_to_stat_entry[v] = m_mask_val_to_stat_entry.size();
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

    // Init stat to 0
    Kokkos::deep_copy(m_stat_field.view_nonconst<Real>(),0.0);

    // We'll use this to decide which of the input field's indices
    // is needed to index the stat field
    auto mask_dim_pos = m_field.layout().dim_idx(m_mask_field.layout().names()[0]);

    int midx;
    view_1d_host<const Real> w;
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
          if (m_use_weight) {
            w = m_weight_field.part_view<const Real>(p);
          }

          // Loop over entries of this part
          for (int i=0; i<fpl.dims()[0]; ++i) {
            midx = m_mask_val_to_stat_entry.at(m(offset+i));

            // Update field (weighted) integral (and weight integral)
            if (m_use_weight) {
              sview(midx) += fview(i)*w(i);
            } else {
              sview(midx) += fview(i);
            }
          }

          // Update offset
          offset += fpl.dims()[0];
        }

        // Global reduction of (weighted) integral
        m_comm.all_reduce(sview.data(),sview.size(),MPI_SUM);
        break;
      }
      case 2:
      {
        // Offsets of i/j part-indices into the unpartitioned field
        int offset_i = 0;
        int offset_j = 0;

        // Init stat to 0
        auto sview = m_stat_field.nd_view_nonconst<Real,2>();
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
            // If i index is the masked one, get stat mask index
            if (mask_dim_pos==0)
              midx = m_mask_val_to_stat_entry.at(m(i+offset_i));

            for (int j=0; j<fpl.dims()[1]; ++j) {
              // If j index is the masked one, get stat mask index
              if (mask_dim_pos==1)
                midx = m_mask_val_to_stat_entry.at(m(j+offset_j));

              // Use non-masked index to access stat,
              // and masked index to access the weight
              if (m_use_weight) {
                if (mask_dim_pos==0) {
                  sview(midx,j+offset_j) += w(i)*fview(i,j);
                } else {
                  sview(i+offset_i,midx) += w(j)*fview(i,j);
                }
              } else {
                // Simply add
                if (mask_dim_pos==0) {
                  sview(midx,j+offset_j) += fview(i,j);
                } else {
                  sview(i+offset_i,midx) += fview(i,j);
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
        auto sview = m_stat_field.nd_view_nonconst<Real,3>();
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
            // If i index is the masked one, get stat mask index
            if (mask_dim_pos==0)
              midx = m_mask_val_to_stat_entry.at(m(i+offset_i));

            for (int j=0; j<fpl.dims()[1]; ++j) {
              // If j index is the masked one, get stat mask index
              if (mask_dim_pos==1)
                midx = m_mask_val_to_stat_entry.at(m(j+offset_j));

              for (int k=0; k<fpl.dims()[2]; ++k) {
                // If k index is the masked one, get stat mask index
                if (mask_dim_pos==2)
                  midx = m_mask_val_to_stat_entry.at(m(k+offset_k));

                // Use non-masked index to access stat,
                // and masked index to access the weight
                if (m_use_weight) {
                  if (mask_dim_pos==0) {
                    sview(midx,j+offset_j,k+offset_k) += w(i)*fview(i,j,k);
                  } else if (mask_dim_pos==1) {
                    sview(i+offset_i,midx,k+offset_k) += w(j)*fview(i,j,k);
                  } else {
                    sview(i+offset_i,j+offset_j,midx) += w(k)*fview(i,j,k);
                  }
                } else {
                  // Simply add
                  if (mask_dim_pos==0) {
                    sview(midx,j+offset_j,k+offset_k) += fview(i,j,k);
                  } else if (mask_dim_pos==1) {
                    sview(i+offset_i,midx,k+offset_k) += fview(i,j,k);
                  } else {
                    sview(i+offset_i,j+offset_j,midx) += fview(i,j,k);
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

    if (m_average) {
      auto wint_v = m_weight_integral.view<Real>();
      switch (fl.rank()) {
        case 1:
        {
          auto sview = m_stat_field.nd_view_nonconst<Real,1>();
          for (int i=0; i<sview.extent_int(0); ++i) {
            sview(i) /= wint_v(i);
          }
          break;
        }
        case 2:
        {
          auto sview = m_stat_field.nd_view_nonconst<Real,2>();
          for (int i=0; i<sview.extent_int(0); ++i) {
            for (int j=0; j<sview.extent_int(1); ++j) {
              if (mask_dim_pos==0) {
                sview(i,j) /= wint_v(i);
              } else {
                sview(i,j) /= wint_v(j);
              }
            }
          }
          break;
        }
        case 3:
        {
          auto sview = m_stat_field.nd_view_nonconst<Real,3>();
          for (int i=0; i<sview.extent_int(0); ++i) {
            for (int j=0; j<sview.extent_int(1); ++j) {
              for (int k=0; k<sview.extent_int(2); ++k) {
                if (mask_dim_pos==0) {
                  sview(i,j,k) /= wint_v(i);
                } else if (mask_dim_pos==1) {
                  sview(i,j,k) /= wint_v(j);
                } else {
                  sview(i,j,k) /= wint_v(k);
                }
              }
            }
          }
          break;
        }
        default:
          EKAT_ERROR_MSG ("Error! [FieldMaskedIntegral] Unsupported field rank (" + std::to_string (fl.rank()) + ")\n");
      }
    }
  }

  // The mask field
  Field         m_mask_field;

  // Map every mask value to an index in [0,N), with N=number_of_mask_values
  std::map<int,int>   m_mask_val_to_stat_entry;
  
  // Optionally, we weigh the integrand by a weight field
  bool          m_use_weight;
  bool          m_average;
  Field         m_weight_field;
  Field         m_weight_integral;
};

} // namespace cldera

#endif // CLDERA_FIELD_MASKED_INTEGRAL_HPP
