#include "cldera_field_masked_integral.hpp"
#include "io/cldera_pnetcdf.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera {

FieldMaskedIntegral::
FieldMaskedIntegral (const ekat::Comm& comm,
                    const ekat::ParameterList& pl)
 : FieldStat(comm,pl)
{
  m_use_weight = m_params.isParameter("weight_field");
  m_average = m_params.get<bool>("average",true);
  m_output_mask_field = m_params.get("output_mask_field",false);
}

std::vector<std::string>
FieldMaskedIntegral::
get_aux_fields_names () const
{
  std::vector<std::string> aux_fnames;
  aux_fnames.push_back("col_gids");
  if (not m_output_mask_field and m_use_weight) {
    aux_fnames.push_back(m_params.get<std::string>("weight_field"));
  }
  return aux_fnames;
}

FieldLayout 
FieldMaskedIntegral::
stat_layout (const FieldLayout& fl) const
{
  if (m_output_mask_field) {
    return m_mask_field.layout();
  }

  const auto& masked_dim = m_mask_field.layout().names()[0];
  auto names = fl.names();
  auto dims  = fl.dims();
  auto pos = fl.dim_idx(masked_dim);
  dims[pos] = m_mask_val_to_stat_entry.size();
  names[pos] = "dim" + std::to_string(dims[pos]);
  return FieldLayout(dims,names);
}

void FieldMaskedIntegral::
set_aux_fields_impl (const std::map<std::string,Field>& fields)
{
  // For sure we have the col gids
  const auto& gids = fields.at("col_gids");
  load_mask_field(gids);

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

    // if (m_weight_field.nparts()!=m_field.nparts()
    if (m_weight_field.nparts()>1) {
      // For simplicity, make w a single-part field
      auto& s = StatFactory::instance();
      ekat::ParameterList pl;
      pl.set("name",m_name + "::set_aux_fields_impl::weight_field");
      auto id = s.create("identity",m_comm,pl);
      id->set_field(m_weight_field);
      id->create_stat_field();
      m_weight_field = id->compute(m_timestamp);
    }
  }

  if (m_average) {
    ekat::ParameterList pl("w_int");
    pl.set("mask_field",m_mask_field.name());
    pl.set("average",false);
    pl.set("mask_file_name",m_params.get<std::string>("mask_file_name"));
    FieldMaskedIntegral w_int_stat(m_comm,pl);
    std::map<std::string,Field> aux_fields;
    aux_fields["mask"] = m_mask_field;
    aux_fields["col_gids"] = gids;
    if (m_use_weight) {
      w_int_stat.set_field(m_weight_field);
    } else {
      m_weight_field = Field("",m_mask_field.layout(),DataAccess::Copy);
      m_weight_field.commit();
      Kokkos::deep_copy(m_weight_field.view_nonconst<Real>(),1);
      w_int_stat.set_field(m_weight_field);
    }
    w_int_stat.set_aux_fields (aux_fields);
    w_int_stat.create_stat_field();
    m_weight_integral = w_int_stat.compute(m_timestamp);
  }

  // First, gather all the mask values we have
  auto data = m_mask_field.data<int>();
  auto size = m_mask_field.layout().size();
  std::set<int> my_mask_vals;
  for (int i=0; i<size; ++i) {
    my_mask_vals.insert(data[i]);
  }

  // Broadcast values, so all procs know all mask values
  std::set<int> mask_vals;
  for (int proc=0; proc<m_comm.size(); ++proc) {
    int nvals = my_mask_vals.size();
    m_comm.broadcast(&nvals,1,proc);
    std::vector<int> this_proc_vals;
    if (proc==m_comm.rank()) {
      for (auto v : my_mask_vals) {
        this_proc_vals.push_back(v);
      }
    } else {
      this_proc_vals.resize(nvals);
    }

    m_comm.broadcast(this_proc_vals.data(),nvals,proc);

    for (auto v : this_proc_vals) {
      mask_vals.insert(v);
    }
  }

  // Then, map each mask value to an index in 0,...,num_mask_values-1
  for (auto v : mask_vals) {
    m_mask_val_to_stat_entry[v] = m_mask_val_to_stat_entry.size();
  }
}

void FieldMaskedIntegral::
compute_impl () {
  if (m_output_mask_field) {
    static bool already_computed = false;
    if (not already_computed) {
      auto s = m_stat_field.view_nonconst<Real>();
      auto m = m_mask_field.view<int>();
      for (int i=0; i<m.extent_int(0); ++i) {
        s[i] = m[i];
      }
      already_computed = true;
    }
    return;
  }
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
void FieldMaskedIntegral::
do_compute_impl ()
{
  const auto& fl = m_field.layout();
  auto m = m_mask_field.view<int>();

  // Init stat to 0
  Kokkos::deep_copy(m_stat_field.view_nonconst<Real>(),0.0);

  // We'll use this to decide which of the input field's indices
  // is needed to index the stat field
  auto mask_dim_pos = m_field.layout().dim_idx(m_mask_field.layout().names()[0]);

  int midx;
  view_1d_host<const Real> w_view;
  if (m_use_weight) {
    w_view = m_weight_field.view<const Real>();
  }
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
          midx = m_mask_val_to_stat_entry.at(m(offset+i));

          // Update field (weighted) integral (and weight integral)
          if (m_use_weight) {
            sview(midx) += fview(i)*w_view(i+offset);
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
                sview(midx,j+offset_j) += fview(i,j)*w_view(i+offset_i);
              } else {
                sview(i+offset_i,midx) += fview(i,j)*w_view(j+offset_j);
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
                  sview(midx,j+offset_j,k+offset_k) += fview(i,j,k)*w_view(i+offset_i);
                } else if (mask_dim_pos==1) {
                  sview(i+offset_i,midx,k+offset_k) += fview(i,j,k)*w_view(j+offset_j);
                } else {
                  sview(i+offset_i,j+offset_j,midx) += fview(i,j,k)*w_view(k+offset_k);
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

void FieldMaskedIntegral::
load_mask_field (const Field& my_col_gids)
{
  const auto& filename = m_params.get<std::string>("mask_file_name");
  auto file = io::pnetcdf::open_file (filename,m_comm,io::pnetcdf::IOMode::Read);
  const auto& mask_name = m_params.get<std::string>("mask_field","mask");
  
  EKAT_REQUIRE_MSG (file->vars.find(mask_name)!=file->vars.end(),
      "Error! Mask field not found in the NC file.\n"
      " - file name: " + filename + "\n"
      " - mask name: " + mask_name + "\n");

  // Ensure gids field has just one part
  Field gids_1p;
  if (my_col_gids.nparts()==1) {
    gids_1p = my_col_gids;
  } else {
    auto& s = StatFactory::instance();
    ekat::ParameterList pl;
    pl.set("name",m_name + "::load_mask_field::col_gids");
    auto id = s.create("identity",m_comm,pl);
    id->set_field(my_col_gids);
    id->create_stat_field();
    gids_1p = id->compute(m_timestamp);
  }

  // Compute min gid (so we get offsets right)
  int min_gid;
  {
    auto& s = StatFactory::instance();
    ekat::ParameterList pl;
    pl.set("name",m_name + "::load_mask_field::min_gid");
    auto min_gid_stat = s.create("global_min",m_comm,pl);
    min_gid_stat->set_field(gids_1p);
    min_gid_stat->create_stat_field();
    min_gid = min_gid_stat->compute(m_timestamp).data<int>()[0];
  }

  // Create offsets
  std::vector<int> offsets;
  const int num_gids = gids_1p.layout().size();
  offsets.reserve(num_gids);
  auto gids = gids_1p.data<int>();
  for (int i=0; i<num_gids; ++i) {
    offsets.push_back(gids[i]-min_gid);
  }

  // Add decomp for pnetcdf
  io::pnetcdf::add_decomp (*file,"ncol",offsets);

  // Check mask var specs
  auto print_dims = [] (const std::vector<std::shared_ptr<const io::pnetcdf::NCDim>>& dims) {
    std::string s;
    if (dims.size()>0) {
      s += dims[0]->name;
      for (size_t i=1; i<dims.size(); ++i) {
        s += dims[i]->name;
      }
    }
    return s;
  };
  const auto& mask_var = file->vars.at(mask_name);
  EKAT_REQUIRE_MSG (mask_var->dims.size()==1,
      "Error! FieldMaskedIntegral requires a 1-dim mask field.\n"
      " - stat name: " + name() + "\n"
      " - mask file: " + filename + "\n"
      " - mask name: " + mask_name + "\n"
      " - mask layout: " + print_dims(mask_var->dims) + "\n");
  const auto& mask_dim = mask_var->dims[0]->name;

  EKAT_REQUIRE_MSG (m_field.layout().has_dim_name(mask_dim),
      "Error! Input field does not have mask field dimension in its layout.\n"
      " - stat name: " + name() + "\n"
      " - mask dim name: " + mask_dim + "\n"
      " - field layout : " + ekat::join(m_field.layout().names(),",") + "\n");

  // Read mask
  m_mask_field = Field(mask_name,gids_1p.layout(),DataAccess::Copy,DataType::IntType);
  m_mask_field.commit();
  io::pnetcdf::read_var(*file,mask_name,m_mask_field.data_nonconst<int>());

  io::pnetcdf::close_file(*file);
}

} // namespace cldera