#include "cldera_field_masked_write.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"
#include "io/cldera_pnetcdf.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera {

FieldMaskedWrite::
FieldMaskedWrite (const ekat::Comm& comm,
                  const ekat::ParameterList& pl)
 : FieldStat(comm,pl)
 , m_write_values(m_params.get<std::vector<Real>>("write_values"))
 , m_write_heights(m_params.get<std::vector<Real>>("write_heights"))
 , m_default_write(m_params.get<Real>("default_write"))
 , m_mask_field_name(m_params.get<std::string>("mask_field"))
{
  m_output_mask_field = m_params.get("output_mask_field",false); // GH: is this needed?

  std::cout << "[masked_write] variables:\n"
            << "  m_mask_field_name = " << m_mask_field_name << "\n"
            << "  m_default_write   = " << m_default_write << "\n"
            << "  m_write_values    = [" << m_write_values << "]\n"
            << "  m_write_heights   = [" << m_write_heights << "]\n";
}

std::vector<std::string>
FieldMaskedWrite::
get_aux_fields_names () const
{
  std::vector<std::string> aux_fnames;
  // we need column GIDs to process data on the native grid compared to the mask file
  aux_fnames.push_back("col_gids");
  // we need the name of the mask field in the corresponding .nc file
  aux_fnames.push_back(m_mask_field_name);
  // we need to read the heights of each pressure level in order to compute the height to write at
  aux_fnames.push_back("zi");
  // we need to read the areas of each column in order to compute the volume to write over
  aux_fnames.push_back("area");

  return aux_fnames;
}

FieldLayout 
FieldMaskedWrite::
stat_layout (const FieldLayout& fl) const
{
  // TODO: this should be a (cols,levs) field
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

void FieldMaskedWrite::
set_aux_fields_impl ()
{
  const auto& gids = m_aux_fields.at("col_gids");
  m_height_field = m_aux_fields.at("zi");
  m_area_field = m_aux_fields.at("area");

  if (m_aux_fields.count(m_mask_field_name)>0) {
    m_mask_field = m_aux_fields.at(m_mask_field_name);
  } else {
    load_mask_field(gids);
  }

  // Setup 0: mask preprocessing
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

  // GH: for write, we assume mask values should already be 0,...,n in order
  // Then, map each mask value to an index in 0,...,num_mask_values-1
  for (auto v : mask_vals) {
    m_mask_val_to_stat_entry[v] = m_mask_val_to_stat_entry.size();
  }

  // Setup 1: identify mask columns


  // Setup 2: identify write heights (separate function)


  // Setup 3: compute write volumes (separate function)
}

void FieldMaskedWrite::
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

  const auto dt   = m_field.data_type();
  const int  rank = m_field.layout().rank();
  std::cout << "[masked_write] field rank = " << rank << std::endl;
  std::cout << "[masked_write] zi rank = " << m_height_field.layout().rank() << std::endl;
  if (dt==DataType::RealType) {
    switch (rank) {
      case 1: return do_compute_impl<Real,1>();
      case 2: return do_compute_impl<Real,2>();
      case 3: return do_compute_impl<Real,3>();
    }
  } else if (dt==DataType::IntType) {
    switch (rank) {
      case 1: return do_compute_impl<int,1>();
      case 2: return do_compute_impl<int,2>();
      case 3: return do_compute_impl<int,3>();
    }
  } else {
    EKAT_ERROR_MSG ("Error! Unexpected/unsupported field data type.\n"
        " - field name: " + m_field.name() + "\n"
        " - field data type: " + e2str(m_field.data_type()) + "\n"
        " - stat name: " + name () + "\n");
  }
}


/*
Set aux fields. This means we have "zi" and "mask" and "col_gids" fields available now.
Using mask, loop over all the cells and determine the IDs of the columns for the injection sites. For each injection site, save a list of those local ids.
Wipe the field and replace it with the default write value
Using those local ids, loop over the mask injection height array. For each entry in the array, loop over all of the corresponding local columns.

for ipart in parts:
  for index in mask:
    

for iheight,height in enumerate(heights):
  num_cols = len(mask_cols[iheight])
  for icol,col in enumerate(mask_cols[iheight]):
    for iz,z in enumerate(zi[col]):
      if z > height:
        injection_level[iheight][icol] = iz
        break

    so2[col][iz] = convert(values[iheight])
*/

template<typename T, int N>
void FieldMaskedWrite::
do_compute_impl ()
{
  std::cout << "[masked_write] called do_compute_impl with N = " << N << std::endl;
  auto sview = m_stat_field.nd_view_nonconst<Real,N>();
  auto mview = m_mask_field.view<int>();

  // Init stat to 0
  Kokkos::deep_copy(sview,0.0);

  const auto& mask_dim_name = m_mask_field.layout().names()[0]; // should be columns
  const int mask_dim = m_field.layout().dim_idx(mask_dim_name); // grab the index for columns
  const int part_dim = m_field.part_dim();

  std::cout << "[masked_write] mask_dim_name = " << mask_dim_name << std::endl;

  // PART 1: construct m_mask_cols
  // Get the columns corresponding to each mask index: e.g. m_mask_cols[i] = [11,12,13] means columns 11,12,13 correspond to write location i
  for (int p=0; p<m_field.nparts(); ++p) {
    auto fpl = m_field.part_layout(p);
    auto fview = m_field.part_nd_view<T,N>(p);

    const int part_offset = m_field.part_offset(p);

    const int mask_dim_offset = mask_dim==part_dim ? part_offset : 0;
    const int mask_dim_ext = fpl.extent(mask_dim_name);
    for (int i=0; i<mask_dim_ext; ++i) {
      auto mval = mview(i+mask_dim_offset);
      auto midx = m_mask_val_to_stat_entry.at(mval);
      std::cout << mval << " ";
    }
    for (int i=0; i<mask_dim_ext; ++i) {
      auto mval = mview(i+mask_dim_offset);
      auto midx = m_mask_val_to_stat_entry.at(mval);
      std::cout << midx << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;


  // NOTE: we operate under the assumption that either
  //  - mask_dim==part_dim: we're integrating over the partitioned dim
  //  - nparts==1: which dim is "partitioned" is arbitrary and irrelevant
  // If we allowed a non-masked dim to be partitioned, then we would have
  // to offset the stat indices also along the non-mask dimensions.
  // That's too many cases, which are likely never needed

  for (int p=0; p<m_field.nparts(); ++p) {
    auto fpl = m_field.part_layout(p);
    auto fview = m_field.part_nd_view<T,N>(p);

    // Offset of this part into the unpartitioned dimension
    const int part_offset = m_field.part_offset(p);

    const int mask_dim_offset = mask_dim==part_dim ? part_offset : 0;
    const int mask_dim_ext = fpl.extent(mask_dim_name);
    for (int i=0; i<mask_dim_ext; ++i) {
      auto mval = mview(i+mask_dim_offset);
      auto midx = m_mask_val_to_stat_entry.at(mval);
      auto w = 1;
      if constexpr (N==1) {
        sview(midx) += fview(i) * w;
      } else {
        auto f_slice = slice(fview,mask_dim,i);
        auto s_slice = slice(sview,mask_dim,midx);
        for (int j=0; j<f_slice.extent_int(0); ++j) {
          if constexpr (N==2) {
            s_slice(j) += f_slice(j) * w;
          } else {
             for (int k=0; k<f_slice.extent_int(1); ++k) {
               s_slice(j,k) += f_slice(j,k) * w;
             }
          }
        }
      }
    }
  }
  track_mpi_all_reduce(m_comm,sview.data(),sview.size(),MPI_SUM,name());
}

void FieldMaskedWrite::
load_mask_field (const Field& my_col_gids)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer (name()+"_load_mask");
  const auto& filename = m_params.get<std::string>("mask_file_name");
  auto file = io::pnetcdf::open_file (filename,m_comm,io::pnetcdf::IOMode::Read);
  
  const auto& mask_name = m_params.get<std::string>("mask_field");
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
      "Error! FieldMaskedWrite requires a 1-dim mask field.\n"
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
  Field mask_field(mask_name,gids_1p.layout(),DataAccess::Copy,DataType::IntType);
  mask_field.commit();
  io::pnetcdf::read_var(*file,mask_name,mask_field.data_nonconst<int>());

  io::pnetcdf::close_file(*file);

  // Store in aux fields. Mark read only, to avoid possible data corruption.
  m_aux_fields[mask_name] = m_mask_field = mask_field.read_only();
  ts.stop_timer (name()+"_load_mask");
}

} // namespace cldera
