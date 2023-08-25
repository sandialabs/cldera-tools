#include "io/cldera_pnetcdf.hpp"
#include "timing/cldera_timing_session.hpp"

#include "cldera_config.h"

#include <ekat/ekat_assert.hpp>
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <pnetcdf.h>

namespace cldera {
namespace io {
namespace pnetcdf {

int UNLIMITED = NC_UNLIMITED;

template<typename T>
void pack_1d (const T* input, const int* dims, int idim, int slice, T* output) {
  EKAT_REQUIRE_MSG (idim==0, "Invalid idim for pack_1d.\n");
  EKAT_REQUIRE_MSG (slice>=0 && slice<dims[0], "Invalid slice for pack_1d.\n");

  *output = input[slice];
}

template<typename T>
void pack_2d (const T* input, const int* dims, int idim, int slice, T* output) {
  EKAT_REQUIRE_MSG (idim==0 || idim==1, "Invalid idim for pack_2d.\n");
  EKAT_REQUIRE_MSG (slice>=0 && slice<dims[idim], "Invalid slice for pack_2d.\n");

  std::vector<int> beg(2,0);
  std::vector<int> end = {dims[0], dims[1]};
  beg[idim] = slice;
  end[idim] = slice+1;
  for (int i=beg[0]; i<end[0]; ++i) {
    for (int j=beg[1]; j<end[1]; ++j,++output) {
      *output = input[i*dims[1]+j];
    }
  }
}

template<typename T>
void pack_3d (const T* input, const int* dims, int idim, int slice, T* output) {
  EKAT_REQUIRE_MSG (idim==0 || idim==1 || idim==2, "Invalid idim for pack_3d.\n");
  EKAT_REQUIRE_MSG (slice>=0 && slice<dims[idim], "Invalid slice for pack_3d.\n");

  std::vector<int> beg(3,0);
  std::vector<int> end = {dims[0], dims[1], dims[2]};
  beg[idim] = slice;
  end[idim] = slice+1;
  for (int i=beg[0]; i<end[0]; ++i) {
    for (int j=beg[1]; j<end[1]; ++j) {
      for (int k=beg[2]; k<end[2]; ++k,++output) {
        *output = input[i*dims[1]*dims[2]+j*dims[2]+k];
      }
    }
  }
}

template<typename T>
void unpack_1d (const T* input, const int* dims, int idim, int slice, T* output) {
  EKAT_REQUIRE_MSG (idim==0, "Invalid idim for unpack_1d.\n");
  EKAT_REQUIRE_MSG (slice>=0 && slice<dims[0], "Invalid slice for unpack_1d.\n");

  output[slice] = *input;
}

template<typename T>
void unpack_2d (const T* input, const int* dims, int idim, int slice, T* output) {
  EKAT_REQUIRE_MSG (idim==0 || idim==1, "Invalid idim for unpack_2d.\n");
  EKAT_REQUIRE_MSG (slice>=0 && slice<dims[idim], "Invalid slice for unpack_2d.\n");

  std::vector<int> beg(2,0);
  std::vector<int> end = {dims[0], dims[1]};
  beg[idim] = slice;
  end[idim] = slice+1;
  for (int i=beg[0]; i<end[0]; ++i) {
    for (int j=beg[1]; j<end[1]; ++j,++input) {
      output[i*dims[1]+j] = *input;
    }
  }
}

template<typename T>
void unpack_3d (const T* input, const int* dims, int idim, int slice, T* output) {
  EKAT_REQUIRE_MSG (idim==0 || idim==1 || idim==2, "Invalid idim for unpack_3d.\n");
  EKAT_REQUIRE_MSG (slice>=0 && slice<dims[idim], "Invalid slice for unpack_3d.\n");

  std::vector<int> beg(3,0);
  std::vector<int> end = {dims[0], dims[1], dims[2]};
  beg[idim] = slice;
  end[idim] = slice+1;
  for (int i=beg[0]; i<end[0]; ++i) {
    for (int j=beg[1]; j<end[1]; ++j) {
      for (int k=beg[2]; k<end[2]; ++k,++input) {
        output[i*dims[1]*dims[2]+j*dims[2]+k] = *input;
      }
    }
  }
}

// If you add support for other types, add specializations of these functions
template<typename T>
nc_type get_nc_type ();
template<typename T>
MPI_Datatype get_io_mpi_dtype ();

template<>
std::string get_io_dtype_name<int> () { return "int"; }
template<>
std::string get_io_dtype_name<long long> () { return "long long"; }
template<>
std::string get_io_dtype_name<float> () { return "float"; }
template<>
std::string get_io_dtype_name<double> () { return "double"; }

template<>
nc_type get_nc_type<int> () { return NC_INT; }
template<>
nc_type get_nc_type<long long> () { return NC_INT64; }
template<>
nc_type get_nc_type<float> () { return NC_FLOAT; }
template<>
nc_type get_nc_type<double> () { return NC_DOUBLE; }
template<>
nc_type get_nc_type<char> () { return NC_CHAR; }

template<>
MPI_Datatype get_io_mpi_dtype<int> () { return MPI_INT; }
template<>
MPI_Datatype get_io_mpi_dtype<long long> () { return MPI_LONG_LONG_INT; }
template<>
MPI_Datatype get_io_mpi_dtype<float> () { return MPI_FLOAT; }
template<>
MPI_Datatype get_io_mpi_dtype<double> () { return MPI_DOUBLE; }

int get_nc_type (const std::string& dtype) {
  if (dtype=="int") {
    return NC_INT;
  } else if (dtype=="long long") {
    return NC_INT64;
  } else if (dtype=="float") {
    return NC_FLOAT;
  } else if (dtype=="double") {
    return NC_DOUBLE;
  } else if (dtype=="char") {
    return NC_CHAR;
  }

  // Not a type
  return NC_NAT;
}

// Helper fcn, for debug/print purposes
std::string dims_str (const std::vector<std::shared_ptr<const NCDim>>& dims) {
  std::vector<std::string> names;
  for (const auto& d : dims) {
    names.push_back(d->name);
  }

  return ekat::join(names,",");
}

std::string dims_str (const std::map<std::string,std::shared_ptr<NCDim>>& dims) {
  std::vector<std::string> names;
  for (const auto& d : dims) {
    names.push_back(d.second->name);
  }
  return ekat::join(names,",");
}

// ==================== QUERY OPS =================== //

bool has_dim (const NCFile& file, const std::string& dname)
{
  return file.dims.find(dname)!=file.dims.end();
}

bool has_var (const NCFile& file,const std::string& vname)
{
  return file.vars.find(vname)!=file.vars.end();
}

bool has_dim (const NCVar& var, const std::string& dname)
{
  for (auto d : var.dims) {
    if (d->name==dname) {
      return true;
    }
  }
  return false;
}


// ==================== FILE OPS ==================== //

std::shared_ptr<NCFile>
open_file (const std::string& fname, const ekat::Comm& comm, const IOMode mode)
{
  int ret;
  auto file = std::make_shared<NCFile>();
  if (mode==IOMode::Write) {
    int ncmode = NC_CLOBBER | NC_64BIT_DATA;
    ret = ncmpi_create(comm.mpi_comm(), fname.c_str(), ncmode, MPI_INFO_NULL, &file->ncid);
  } else {
    int ncmode = mode==IOMode::Append ? NC_WRITE : NC_NOWRITE;
    ret = ncmpi_open(comm.mpi_comm(), fname.c_str(), ncmode, MPI_INFO_NULL, &file->ncid);
  }

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not open file.\n"
      "   - file name: " + fname + "\n"
      "   - open mode: " + e2str(mode) + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file->name = fname;
  file->mode = mode;
  file->comm = comm;

  if (mode==IOMode::Read or mode==IOMode::Append) {
    // Populate file with dims and vars
    char name [NC_MAX_NAME];
    int ndims, nvars, ngatts, unlimited;
    ret = ncmpi_inq(file->ncid, &ndims, &nvars, &ngatts, &unlimited);
    EKAT_REQUIRE_MSG (ret==NC_NOERR,
        "Error! Could not inquire counters.\n"
        "   - file name: " + fname + "\n"
        "   - err code : " + std::to_string(ret) + "\n");

    // Read dims
    MPI_Offset dimlen;
    std::shared_ptr<NCDim> dim;
    std::map<int,std::string> dim_id2name;
    for (int idim=0; idim<ndims; ++idim) {
      dim = std::make_shared<NCDim>();
      ret = ncmpi_inq_dim(file->ncid, idim, name, &dimlen);
      EKAT_REQUIRE_MSG (ret==NC_NOERR,
          "Error! Could not inquire dimension.\n"
          "   - file name: " + fname + "\n"
          "   - dim id   : " + std::to_string(idim) + "\n"
          "   - err code : " + std::to_string(ret) + "\n");
      dim->ncid = idim;
      dim->name = name;
      dim->len = dim->glen = dimlen;
      file->dims.emplace(name,dim);
      dim_id2name[idim] = name;
    }

    // Read vars
    std::shared_ptr<NCVar> var;
    nc_type dtype;
    int var_ndims,var_natts;
    std::vector<int> var_dimids;
    for (int ivar=0; ivar<nvars; ++ivar) {
      var = std::make_shared<NCVar>();

      ret = ncmpi_inq_varndims(file->ncid, ivar, &var_ndims);
      EKAT_REQUIRE_MSG (ret==NC_NOERR,
          "Error! Could not inquire variable's number of dims.\n"
          "   - file name: " + fname + "\n"
          "   - var id   : " + std::to_string(ivar) + "\n"
          "   - err code : " + std::to_string(ret) + "\n");
      var_dimids.resize(var_ndims);

      ret = ncmpi_inq_var(file->ncid, ivar, name, &dtype, &var_ndims,
                          var_dimids.data(), &var_natts);
      EKAT_REQUIRE_MSG (ret==NC_NOERR,
          "Error! Could not inquire variable.\n"
          "   - file name: " + fname + "\n"
          "   - var id   : " + std::to_string(ivar) + "\n"
          "   - err code : " + std::to_string(ret) + "\n");

      var->name = name;
      var->ncid = ivar;
      if (dtype==NC_DOUBLE) {
        var->dtype = "double";
      } else if (dtype==NC_FLOAT) {
        var->dtype = "float";
      } else if (dtype==NC_INT) {
        var->dtype = "int";
      } else if (dtype==NC_INT64) {
        var->dtype = "long long";
      } else {
        EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type: " + std::to_string(dtype) + "\n");
      }
      var->size = 1;
      for (int idim=0; idim<var_ndims; ++idim) {
        auto dname = dim_id2name.at(var_dimids[idim]);
        var->dims.push_back(file->dims.at(dname));
        var->size *= var->dims.back()->len;
      }

      if (var->dims[0]->name=="time") {
        var->nrecords = file->dims.at("time")->len;
      }

      file->vars.emplace(name,var);
    }
    file->enddef = true;
  } else {
    // Add time unlimited dim
    add_dim (*file,"time",NC_UNLIMITED);
  }

  return file;
}

void close_file (NCFile& file)
{
  int ret = ncmpi_close(file.ncid);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not close NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file.ncid = -1;
  file.mode = IOMode::Invalid;
  file.name = "";
  file.dims.clear();
  file.vars.clear();
}

void enddef (NCFile& file)
{
  int ret = ncmpi_enddef(file.ncid);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not end define mode on NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file.enddef = true;
}
void redef (NCFile& file)
{
  int ret = ncmpi_redef(file.ncid);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not re-enter define mode on NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file.enddef = false;
}

void add_dim (NCFile& file, const std::string& dname, const int len, const bool partitioned)
{
  EKAT_REQUIRE_MSG (not file.enddef,
      "Error! Cannot add dimensions after enddef was called.\n"
      "  - file name: " + file.name + "\n"
      "  - dim name : " + dname + "\n"
      "  - dim len  : " + std::to_string(len) + "\n");

  if (has_dim(file,dname)) {
    // Check that the dimension is being re-registered in the same way
    const auto& dim = file.dims.at(dname);
    auto mpi_sum = [&] (const int val) -> int {
      int gval;
      file.comm.all_reduce(&val,&gval,1,MPI_SUM);
      return gval;
    };
    EKAT_REQUIRE_MSG ( partitioned==dim->is_partitioned,
        "Error! Could not add dimension to NC file. Dimension already added with different partitioned flag.\n"
        "   - file name: " + file.name + "\n"
        "   - dim name : " + dname + "\n"
        "   - dim partitioned : " + (dim->is_partitioned ? "yes" : "no") + "\n"
        "   - input partitioned : " + (partitioned ? "yes" : "no") + "\n");

    EKAT_REQUIRE_MSG ( dim->len==len,
        "Error! Could not add dimension to NC file. Dimension already added with different local length.\n"
        "   - file name: " + file.name + "\n"
        "   - dim name : " + dname + "\n"
        "   - dim len  : " + std::to_string(dim->len) + "\n"
        "   - input len: " + std::to_string(len) + "\n");
    EKAT_REQUIRE_MSG ( not partitioned || dim->glen==mpi_sum(len),
        "Error! Could not add partitioned dimension to NC file. Dimension already added with different global length.\n"
        "   - file name : " + file.name + "\n"
        "   - dim name  : " + dname + "\n"
        "   - dim glen  : " + std::to_string(dim->glen) + "\n"
        "   - input glen: " + std::to_string(mpi_sum(len)) + "\n");
    return;
  }

  auto dim = std::make_shared<NCDim>();
  dim->name = dname;
  dim->len = dim->glen = len;
  if (partitioned) {
    file.comm.all_reduce(&dim->glen,1,MPI_SUM);
    dim->is_partitioned = true;
  }
  int ret = ncmpi_def_dim(file.ncid, dname.c_str(),dim->glen,&dim->ncid);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not add dimension to NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - dim name : " + dname + "\n"
      "   - dim len  : " + std::to_string(dim->len) + "\n"
      "   - dim glen : " + std::to_string(dim->glen) + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file.dims[dname] = dim;
}

void add_decomp (      NCFile& file, 
                 const std::string& dim_name,
                 const std::vector<int>& entries)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (file.enddef,
      "Error! Cannot add decomposition until enddef has been called.\n"
      "  - file name: " + file.name + "\n"
      "  - dim name : " + dim_name + "\n");

  EKAT_REQUIRE_MSG (file.dims.find(dim_name)!=file.dims.end(),
      "Error! Invalid decomp dimension.\n"
      "   - file name  : " + file.name + "\n"
      "   - file dims  : (" + dims_str(file.dims) + ")\n"
      "   - decomp dim : " + dim_name + "\n");
  auto dim = file.dims.at(dim_name);

  EKAT_REQUIRE_MSG (not dim->decomp_set || dim->entries==entries,
      "Error! Decomposition was already added for this dimension, with different entries.\n"
      "   - file name : " + file.name + "\n"
      "   - dim name  : " + dim_name + "\n"
      "   - decomp entries: [" + ekat::join(dim->entries," ") + "]\n"
      "   - input entries: [" + ekat::join(entries," ") + "]\n");

  if (dim->decomp_set) {
    // We already set the decomp for this dimension
    return;
  }

  // Count entries, so we can check correctness, and, for Read mode,
  // we can adjust the local length of dimensions
  const auto& comm = file.comm;
  int count = entries.size();
  int gcount;
  comm.all_reduce(&count,&gcount,1,MPI_SUM);

  EKAT_REQUIRE_MSG (gcount==dim->glen,
      "Error! Invalid global count for decomposition.\n"
      "   - file name    : " + file.name + "\n"
      "   - decomp dim   : " + dim_name + "\n"
      "   - global count : " + std::to_string(gcount) + "\n"
      "   - dim glength  : " + std::to_string(dim->glen) + "\n");

  // When reading from file, all dims are created as not partitioned.
  // The user is allowed to tell the IO layer that a certain dim is
  // partitioned in a subsequent phase. On the other hand, when we
  // open a file for writing, the dimension must have been registered
  // as partitioned.
  auto read = file.mode==IOMode::Read or file.mode==IOMode::Append;
  if (read) {
    dim->is_partitioned = true;

    // Update the local extent of the dimension, and recompute
    // the variables sizes/dimlens.
    dim->len = count;
    for (auto v : file.vars) {
      v.second->compute_extents();
    }
  } else {
    EKAT_REQUIRE_MSG (dim->is_partitioned,
        "Error! This dimension was not marked as partitioned.\n"
        "   - file name : " + file.name + "\n"
        "   - dim name  : " + dim_name + "\n");

    EKAT_REQUIRE_MSG (count==dim->len,
        "Error! Invalid local count for decomposition.\n"
        "   - file name    : " + file.name + "\n"
        "   - decomp dim   : " + dim_name + "\n"
        "   - mpi rank     : " + std::to_string(comm.rank()) + "\n"
        "   - local count  : " + std::to_string(count) + "\n"
        "   - dim length   : " + std::to_string(dim->len) + "\n");
  }

  // Ensure indices are 1-based (for netcdf)
  int min = *std::min_element(entries.begin(),entries.end());
  comm.all_reduce(&min,1,MPI_MIN);
  EKAT_REQUIRE_MSG (min==0,
      "Error! Invalid index base for dim decomp enries.\n"
      "   - file name  : " + file.name + "\n"
      "   - decomp dim : " + dim_name + "\n"
      "   - min entry  : " + std::to_string(min) + "\n");

  dim->decomp_set = true;
  dim->entries = entries;

  // Loop over vars, to add a decomp for each layout
  for (auto it : file.vars) {
    auto var = it.second;
    auto dims = var->dims;
    int rank = dims.size();
    
    // Check if this var contains the decomposed dim,
    // and keep track of the decomposed dim index
    int idecomp = -1;
    for (int i=0; i<rank; ++i) {
      if (dims[i]->name==dim_name) {
        idecomp=i;
        break;
    }}
          
    if (idecomp>=0) {
      // Create decomp name (and check that there is only 1 decomposed dim)
      std::string name;
      int ndecomps = 0;
      for (auto d : dims) {
        name += d->name + "-";
        if (d->is_partitioned) {
          ++ndecomps;
        }
      }
      name.pop_back();

      EKAT_REQUIRE_MSG (ndecomps==1,
          "Error! A variable cannot be decomposed along more than one dimension!\n");

      if (file.decomps.find(name)==file.decomps.end()) {
        // This layout decomposition is new, add it in the file
        auto& decomp = file.decomps[name];

        // Create a decomp object, and fill its data
        decomp = std::make_shared<IODecomp>();
        decomp->name = name;
        decomp->layout = dims;
        decomp->dim = file.dims.at(dim_name);
        decomp->dim_idx = idecomp;
        decomp->hyperslab_size = var->size / decomp->dim->len;

        // Make buf large enough for any data type
        decomp->buf.resize(sizeof(double)*decomp->hyperslab_size*sizeof(int));

        file.decomps[name] = decomp;
      }
      auto decomp = file.decomps.at(name);
      var->decomp = decomp;
    }
  }
}

void add_var (      NCFile& file,
              const std::string& vname,
              const std::string& dtype,
                    std::vector<std::string> dims,
              const bool time_dep)
{
  EKAT_REQUIRE_MSG (not file.enddef,
      "Error! Cannot add variables after enddef was called.\n"
      "  - file name: " + file.name + "\n"
      "  - var name : " + vname + "\n"
      "  - var dims : (" + ekat::join(dims,",") + ")\n");

  if (has_var(file,vname)) {
    // Check that they're the same, then early return
    const auto& var = *file.vars.at(vname);
    EKAT_REQUIRE_MSG (dtype==var.dtype,
        "Error! Could not add variable to NC file. Variable already added with a different data type.\n"
        "   - file name: " + file.name + "\n"
        "   - var name : " + vname + "\n"
        "   - var dtype : " + var.dtype + "\n"
        "   - input dtype : " + dtype + "\n");
    EKAT_REQUIRE_MSG ( ekat::join(dims,",") == dims_str(var.dims),
        "Error! Could not add variable to NC file. Variable already added with different dimensions.\n"
        "   - file name: " + file.name + "\n"
        "   - var name : " + vname + "\n"
        "   - var dims : (" + dims_str(file.vars.at(vname)->dims) + ")\n"
        "   - input dims : (" + ekat::join(dims,",") + ")\n");

    EKAT_REQUIRE_MSG (time_dep==var.has_time(),
        "Error! Could not add variable to NC file. Variable already added with a different time dependency.\n"
        "   - file name: " + file.name + "\n"
        "   - var name : " + vname + "\n"
        "   - var time dep : " + (var.has_time() ? "yes" : "no") + "\n"
        "   - input time dep : " + (time_dep ? "yes" : "no") + "\n");

    return;
  }

  auto it_time = std::find(dims.begin(),dims.end(),"time");
  EKAT_REQUIRE_MSG (it_time==dims.end() || it_time==dims.begin(),
      "Error! Time dimension must be the first!\n"
      "   - file name: " + file.name + "\n"
      "   - var name : " + vname + "\n"
      "   - var dims : (" + dims_str(file.vars.at(vname)->dims) + ")\n");

  EKAT_REQUIRE_MSG (time_dep || it_time==dims.end(),
      "Error! Var was not supposed to be time dependent, but time is one of its dimensions.\n"
      "   - file name: " + file.name + "\n"
      "   - var name : " + vname + "\n"
      "   - var dims : (" + dims_str(file.vars.at(vname)->dims) + ")\n");

  // Add time dimension if needed
  if (time_dep && it_time==dims.end()) {
    dims.insert(dims.begin(),"time");
  }
  EKAT_REQUIRE_MSG (dims.size()>0,
      "Error! Cannot add a non-time dependent variable that has no dimensions.\n"
        "   - file name : " + file.name + "\n"
        "   - var name  : " + vname + "\n");

  auto var = std::make_shared<NCVar>();
  var->name = vname;
  var->dtype = dtype;
  var->nrecords = 0;
  std::vector<int> dims_ids;

  var->size = 1;
  for (const auto& d : dims) {
    EKAT_REQUIRE_MSG (has_dim(file,d),
        "Error! Could not add variable to NC file. Dimension '" + d + "' not in file.\n"
        "   - file name : " + file.name + "\n"
        "   - var name  : " + vname + "\n"
        "   - var_dims  : (" + ekat::join(dims,",") + ")\n"
        "   - time dep  : " + (time_dep ? "yes" : "no") + "\n"
        "   - file dims : (" + dims_str(file.dims) + ")\n");
    dims_ids.push_back(file.dims.at(d)->ncid);

    auto dim = file.dims.at(d);
    var->dims.push_back(file.dims.at(d));
    var->size *= dim->len;
  }

  var->compute_extents ();

  int nc_dtype;
  if (dtype=="double") {
    nc_dtype = NC_DOUBLE;
  } else if (dtype=="float") {
    nc_dtype = NC_FLOAT;
  } else if (dtype=="int") {
    nc_dtype = NC_INT;
  } else if (dtype=="long long") {
    nc_dtype = NC_INT64;
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type: " + dtype + "\n");
  }

  int ret = ncmpi_def_var(file.ncid,vname.c_str(),nc_dtype,dims.size(),dims_ids.data(),&var->ncid);
  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not add variable to NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - var name : " + vname + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file.vars[vname] = var;
}

void NCVar::compute_extents () {
  size = 1;
  dimlens.reserve(has_time() ? dims.size()-1 : dims.size());
  for (const auto& dim : dims) {
    if (dim->name!="time") {
      size *= dim->len;
      dimlens.push_back(dim->len);
    }
  }
}

void add_time (      NCFile& file,
               const std::string& dtype)
{
  add_dim(file,"time",UNLIMITED);
  add_var(file,"time",dtype,{"time"},true);
}

// ================== READ OPS ================ //

template<typename T>
void read_var (const NCFile& file, const std::string& vname,
                     T* const  data, const int record)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("io::read_var");

  EKAT_REQUIRE_MSG (file.vars.find(vname)!=file.vars.end(),
      "Error! Variable not found in output NC file.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + vname + "\n");

  auto var = file.vars.at(vname);
  EKAT_REQUIRE_MSG (var->dtype==get_io_dtype_name<T>(),
      "Error! Incorrect data type.\n"
      "  - file name : " + file.name + "\n"
      "  - var name   : " + var->name + "\n"
      "  - var dtype  : " + var->dtype + "\n"
      "  - input type : " + get_io_dtype_name<T>() + "\n");

  EKAT_REQUIRE_MSG (data!=nullptr,
      "Error! Invalid data pointer.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + vname + "\n");

  const auto& dims = var->dims;
  std::vector<MPI_Offset> start(dims.size()), count(dims.size());
  const bool has_time = has_dim(*var,"time");
  EKAT_REQUIRE_MSG (record<=0 || has_time,
      "Error! Record specified for non-time dependent variable.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + var->name + "\n"
      "  - var dims  : " + dims_str(dims) + "\n"
      "  - record    : " + std::to_string(record) + "\n");

  EKAT_REQUIRE_MSG (!has_time || record<=var->nrecords,
      "Error! Record out of bounds.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + var->name + "\n"
      "  - var nrecords  : " + std::to_string(var->nrecords) + "\n"
      "  - requested record    : " + std::to_string(record) + "\n");

  int first_non_time_idx = has_time ? 1 : 0;
  if (has_time) {
    start[0] = record>=0 ? record : var->nrecords-1;
    count[0] = 1;
  }
  for (size_t idim=first_non_time_idx; idim<dims.size(); ++idim) {
    start[idim] = 0;
    count[idim] = dims[idim]->len;
  }

  int ret;

  auto mpi_dtype = get_io_mpi_dtype<T>();
  // If partitioned, write one decomposed-dim entry at a time
  if (var->decomp) {
    auto decomp = var->decomp;
    auto dim = decomp->dim;

    // This dim_idx is for pack_Nd, which only uses the non-time dims of the array
    int non_time_dim_offset = var->has_time() ? 1 : 0;
    int dim_idx = decomp->dim_idx-non_time_dim_offset;
    int num_nontime_dims = var->dims.size() - non_time_dim_offset;

    // We read only one entry at a time along the decomp dim
    count[decomp->dim_idx] = 1;

    T* buf = reinterpret_cast<T*>(decomp->buf.data());

    ret = ncmpi_begin_indep_data(file.ncid);
    EKAT_REQUIRE_MSG (ret==NC_NOERR,
        "Error! Could not begin independent data mode for read.\n"
        "  - file name : " + file.name + "\n"
        "  - var name  : " + var->name + "\n"
        "  - err code : " + std::to_string(ret) + "\n");
    for (unsigned i=0; i<dim->entries.size(); ++i) {
      int gid = dim->entries[i];
      start[decomp->dim_idx] = gid;
      ret = ncmpi_get_vara(file.ncid,var->ncid,
                           start.data(),count.data(),
                           buf,decomp->hyperslab_size,mpi_dtype);
      switch (num_nontime_dims) {
        case 1: unpack_1d(buf,var->dimlens.data(),dim_idx,i,data); break;
        case 2: unpack_2d(buf,var->dimlens.data(),dim_idx,i,data); break;
        case 3: unpack_3d(buf,var->dimlens.data(),dim_idx,i,data); break;
        default:
          EKAT_ERROR_MSG ("Error! Max number of non-time dims for decomposed vars is 3.\n"
                          " - file name: " + file.name + "\n"
                          " - var name : " + var->name + "\n"
                          " - var dims : " + dims_str(var->dims) + "\n");
      }
#ifdef CLDERA_DEBUG
      EKAT_REQUIRE_MSG (ret==NC_NOERR,
          "Error! Could not read decomposed variable.\n"
          "  - file name : " + file.name + "\n"
          "  - var name  : " + var->name + "\n"
          "  - err code : " + std::to_string(ret) + "\n");
#endif
    }
    ret = ncmpi_end_indep_data(file.ncid);
    EKAT_REQUIRE_MSG (ret==NC_NOERR,
        "Error! Could not end independent data mode for read.\n"
        "  - file name : " + file.name + "\n"
        "  - var name  : " + var->name + "\n"
        "  - err code : " + std::to_string(ret) + "\n");
  } else {
    ret = ncmpi_get_vara_all(file.ncid,var->ncid,
                             start.data(),count.data(),
                             data,var->size,mpi_dtype);
#ifdef CLDERA_DEBUG
    EKAT_REQUIRE_MSG (ret==NC_NOERR,
        "Error! Could not read non-decomposed variable.\n"
        "  - file name : " + file.name + "\n"
        "  - var name  : " + var->name + "\n"
        "  - err code : " + std::to_string(ret) + "\n");
#endif
  }
  ts.stop_timer("io::read_var");
}

// Instantiations 
template void read_var (const NCFile& file, const std::string& vname,
                        int* const data, const int record);
template void read_var (const NCFile& file, const std::string& vname,
                        long long* const data, const int record);
template void read_var (const NCFile& file, const std::string& vname,
                        float* const data, const int record);
template void read_var (const NCFile& file, const std::string& vname,
                        double* const data, const int record);

// ====================== WRITE OPS ==================== //

template<typename T>
void write_var (const NCFile& file,const std::string& vname,
                const T* const  data)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("io::write_var");
  EKAT_REQUIRE_MSG (file.vars.find(vname)!=file.vars.end(),
      "Error! Variable not found in output NC file.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + vname + "\n");

  auto var = file.vars.at(vname);
  EKAT_REQUIRE_MSG (var->dtype==get_io_dtype_name<T>(),
      "Error! Incorrect data type.\n"
      "  - var name   : " + var->name + "\n"
      "  - var dtype  : " + var->dtype + "\n"
      "  - input type : " + get_io_dtype_name<T>() + "\n");

  EKAT_REQUIRE_MSG (data!=nullptr,
      "Error! Invalid data pointer.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + vname + "\n");

  const auto& dims = var->dims;
  const int ndims = dims.size();
  std::vector<MPI_Offset> start(ndims), count(ndims);
  const int first_non_time_idx = var->has_time() ? 1 : 0;
  const int has_non_time_dims = first_non_time_idx<ndims;
  if (var->has_time()) {
    start[0] = var->nrecords;
    count[0] = 1;
  }
  if (has_non_time_dims) {
    for (size_t idim=first_non_time_idx; idim<dims.size(); ++idim) {
      start[idim] = 0;
      count[idim] = dims[idim]->len;
    }
  }

  int ret;

  ret = ncmpi_begin_indep_data(file.ncid);
  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not begin independent data mode for write.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + var->name + "\n"
      "  - err code : " + std::to_string(ret) + "\n");

  auto mpi_dtype = get_io_mpi_dtype<T>();

  // If a decomposition was provided that impacts this var, we need to do things differently
  if (var->decomp) {
    auto decomp = var->decomp;
    auto dim = decomp->dim;

    // This dim_idx is for pack_Nd, which only uses the non-time dims of the array
    int non_time_dim_offset = var->has_time() ? 1 : 0;
    int dim_idx = decomp->dim_idx-non_time_dim_offset;
    int num_nontime_dims = var->dims.size() - non_time_dim_offset;

    // We write only one entry at a time along the decomp dim
    count[decomp->dim_idx] = 1;

    T* buf = reinterpret_cast<T*>(decomp->buf.data());
    for (unsigned i=0; i<dim->entries.size(); ++i) {
      int gid = dim->entries[i];
      start[decomp->dim_idx] = gid;
      switch (num_nontime_dims) {
        case 1: pack_1d(data,var->dimlens.data(),dim_idx,i,buf); break;
        case 2: pack_2d(data,var->dimlens.data(),dim_idx,i,buf); break;
        case 3: pack_3d(data,var->dimlens.data(),dim_idx,i,buf); break;
        default:
          EKAT_ERROR_MSG ("Error! Max number of non-time dims for decomposed vars is 3.\n"
                          " - file name: " + file.name + "\n"
                          " - var name : " + var->name + "\n"
                          " - var dims : " + dims_str(var->dims) + "\n");
      }
      ret = ncmpi_put_vara(file.ncid,var->ncid,
                           start.data(),count.data(),
                           buf,decomp->hyperslab_size,mpi_dtype);
#ifdef CLDERA_DEBUG
      EKAT_REQUIRE_MSG (ret==NC_NOERR,
          "Error! Could not write partitioned variable.\n"
          "  - file name : " + file.name + "\n"
          "  - var name  : " + var->name + "\n"
          "  - err code : " + std::to_string(ret) + "\n");
#endif
    }
  } else {
    // No partitioned data: simply write to file from root rank.
    // WARNING: if the data is not partitioned, we assume that all
    //          ranks store the same data, and therefore we can
    //          pick any rank to do the write.
    if (file.comm.am_i_root()) {
      ret = ncmpi_put_vara(file.ncid,var->ncid,
                           start.data(),count.data(),
                           data,var->size,mpi_dtype);
#ifdef CLDERA_DEBUG
      EKAT_REQUIRE_MSG (ret==NC_NOERR,
          "Error! Could not write non-partitioned variable.\n"
          "  - file name : " + file.name + "\n"
          "  - var name  : " + var->name + "\n"
          "  - err code : " + std::to_string(ret) + "\n");
#endif
    }
  }

  ret = ncmpi_end_indep_data(file.ncid);
  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not end independent data mode for write.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + var->name + "\n"
      "  - err code : " + std::to_string(ret) + "\n");

  // Update number of records
  ++var->nrecords;
  ts.stop_timer("io::write_var");
}

// Instantiations
template void write_var (const NCFile& file, const std::string& vname,
                         const int* const data);
template void write_var (const NCFile& file, const std::string& vname,
                         const long long* const data);
template void write_var (const NCFile& file, const std::string& vname,
                         const float* const data);
template void write_var (const NCFile& file, const std::string& vname,
                         const double* const data);

// ======================== ATTRIBUTES =========================== //

template<typename T>
void set_att_v (const NCFile& file,
                const std::string& att_name,
                const std::string& var_name,
                const std::vector<T>& data)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("io::write_att");

  int varid;
  if (var_name=="NC_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    EKAT_REQUIRE_MSG (file.vars.find(var_name)!=file.vars.end(),
        "Error! Could not write variable attribute. Variable not found.\n"
        "  - file name : " + file.name + "\n"
        "  - var name  : " + var_name + "\n");

    varid = file.vars.at(var_name)->ncid;
  }

  int ret = ncmpi_put_att (file.ncid,varid,att_name.c_str(),
                           get_nc_type<T>(),data.size(),data.data());

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not write attribute.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + var_name + "\n"
      "  - err code : " + std::to_string(ret) + "\n");

  ts.stop_timer("io::write_att");
}

// Instantiations
template void set_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                         const std::vector<int>&  data);
template void set_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                         const std::vector<long long>& data);
template void set_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                         const std::vector<float>& data);
template void set_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                         const std::vector<double>& data);
template void set_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                         const std::vector<char>& data);

template<typename T>
void get_att_v (const NCFile& file,
                const std::string& att_name,
                const std::string& var_name,  // use "NC_GLOBAL" for global attributes
                      std::vector<T>& data)
{
  auto& ts = timing::TimingSession::instance();
  ts.start_timer("io::read_att");

  int varid;
  if (var_name=="NC_GLOBAL") {
    varid = NC_GLOBAL;
  } else {
    EKAT_REQUIRE_MSG (file.vars.find(var_name)!=file.vars.end(),
        "Error! Could not write variable attribute. Variable not found.\n"
        "  - file name : " + file.name + "\n"
        "  - var name  : " + var_name + "\n");

    varid = file.vars.at(var_name)->ncid;
  }

  int ret;
  MPI_Offset len;
  ret = ncmpi_inq_attlen(file.ncid, varid, att_name.c_str(), &len);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not retrieve attribute length.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + var_name + "\n"
      "  - err code : " + std::to_string(ret) + "\n");

  data.resize(len);
  ret = ncmpi_get_att (file.ncid,varid,att_name.c_str(),data.data());

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not read attribute.\n"
      "  - file name : " + file.name + "\n"
      "  - var name  : " + var_name + "\n"
      "  - err code : " + std::to_string(ret) + "\n");
  ts.stop_timer("io::read_att");
}

// Instantiations
template void get_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                               std::vector<int>& data);
template void get_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                               std::vector<long long>& data);
template void get_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                               std::vector<float>& data);
template void get_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                               std::vector<double>& data);
template void get_att_v (const NCFile& file,
                         const std::string& att_name,
                         const std::string& var_name,
                               std::vector<char>& data);

} // namespace pnetcdf
} // namespace io
} // namespace cldera
