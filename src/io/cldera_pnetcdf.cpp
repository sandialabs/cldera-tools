#include "io/cldera_pnetcdf.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <pnetcdf.h>

namespace cldera {

// ==================== FILE OPS ==================== //

std::shared_ptr<NCFile>
open_file (const std::string& fname, const ekat::Comm& comm, const IOMode mode)
{
  EKAT_REQUIRE_MSG (mode==IOMode::Write || mode==IOMode::Read,
      "Error! Invalid IO mode '" + e2str(mode) + "'.\n");

  int ret;
  auto file = std::make_shared<NCFile>();
  std::string mode_str;
  if (mode==IOMode::Write) {
    int ncmode = NC_CLOBBER | NC_64BIT_DATA;
    mode_str = "write";
    ret = ncmpi_create(comm.mpi_comm(), fname.c_str(), ncmode, MPI_INFO_NULL, &file->ncid);
  } else {
    mode_str = "write";
    int ncmode = NC_NOWRITE;
    ret = ncmpi_open(comm.mpi_comm(), fname.c_str(), ncmode, MPI_INFO_NULL, &file->ncid);
  }

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not open file.\n"
      "   - file name: " + fname + "\n"
      "   - open mode: " + mode_str + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file->name = fname;
  file->mode = mode;
  file->comm = comm;

  if (mode==IOMode::Read) {
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
      dim->len = dimlen;
      file->dims.emplace(name,dim);
      file->dim_id2name[idim] = name;
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
        var->dtype = "int64";
      } else {
        EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type: " + std::to_string(dtype) + "\n");
      }
      for (int idim=0; idim<var_ndims; ++idim) {
        auto dname = file->dim_id2name.at(var_dimids[idim]);
        var->dims.push_back(file->dims.at(dname));
      }

      file->vars.emplace(name,var);
    }
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
  file.decomps.clear();
  file.dim_id2name.clear();
}

void enddef (const NCFile& file)
{
  int ret = ncmpi_enddef(file.ncid);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not end define mode on NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - err code : " + std::to_string(ret) + "\n");
}
void redef   (const NCFile& file)
{
  int ret = ncmpi_redef(file.ncid);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not re-enter define mode on NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - err code : " + std::to_string(ret) + "\n");
}

// =================== DIM OPS =================== //

void add_dim (NCFile& file, const std::string& dname, const int len)
{
  EKAT_REQUIRE_MSG (not has_dim(file,dname),
      "Error! Could add dimension to NC file. Dimension already added.\n"
      "   - file name: " + file.name + "\n"
      "   - dim name : " + dname + "\n");

  auto dim = std::make_shared<NCDim>();
  dim->name = dname;
  dim->len = len;
  int ret = ncmpi_def_dim(file.ncid, dname.c_str(),len,&dim->ncid);

  EKAT_REQUIRE_MSG (ret==NC_NOERR,
      "Error! Could not add dimension to NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - dim name : " + dname + "\n"
      "   - dim len  : " + std::to_string(ret) + "\n"
      "   - err code : " + std::to_string(ret) + "\n");

  file.dims[dname] = dim;
}

bool has_dim (const NCFile& file, const std::string& dname)
{
  return file.dims.find(dname)!=file.dims.end();
}

int get_dim_len (const NCFile& file, const std::string& dname)
{
  EKAT_REQUIRE_MSG (has_dim(file,dname),
      "Error! Dimension not found in NC file.\n"
      "   - file name: " + file.name + "\n"
      "   - dim name : " + dname + "\n");

  return file.dims.at(dname)->len;
}

// ====================== VAR OPS ===================== //

void add_var (      NCFile& file,
              const std::string& vname,
              const std::string& dtype,
              const std::vector<std::string>& dims)
{
  // Useful for error messages
  std::string var_dims_str;
  for (const auto& d : dims) {
    var_dims_str += d;
    var_dims_str += " ";
  }
  std::string file_dims_str;
  for (const auto& d : file.dims) {
    file_dims_str += d.second->name;
    file_dims_str += " ";
  }

  EKAT_REQUIRE_MSG (not has_var(file,vname),
      "Error! Could not add variable to NC file. Variable already added.\n"
      "   - file name: " + file.name + "\n"
      "   - var name : " + vname + "\n"
      "   - var dims : " + var_dims_str + "\n");

  auto var = std::make_shared<NCVar>();
  var->name = vname;
  var->dtype = dtype;
  std::vector<int> dims_ids;
  dims_ids.reserve(dims.size());
  for (const auto& d : dims) {
    EKAT_REQUIRE_MSG (has_dim(file,d),
        "Error! Could add variable to NC file. Variable dimension not in file.\n"
        "   - file name : " + file.name + "\n"
        "   - var name  : " + vname + "\n"
        "   - var dims  : " + var_dims_str + "\n"
        "   - file dims : " + file_dims_str + "\n");
    dims_ids.push_back(file.dims.at(d)->ncid);

    var->dims.push_back(file.dims.at(d));
  }

  int nc_dtype;
  if (dtype=="double") {
    nc_dtype = NC_DOUBLE;
  } else if (dtype=="float") {
    nc_dtype = NC_FLOAT;
  } else if (dtype=="int") {
    nc_dtype = NC_INT;
  } else if (dtype=="int64") {
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

void add_decomp (      NCFile& file, 
                 const std::string& dim_name,
                 const int count)
{
  // Useful for error messages
  std::string file_dims_str;
  for (const auto& d : file.dims) {
    file_dims_str += d.first;
    file_dims_str += " ";
  }

  const auto& comm = file.comm;
  const int rank = comm.rank();

  // Sanity checks
  auto dim_it = file.dims.find(dim_name);
  EKAT_REQUIRE_MSG (dim_it!=file.dims.end(),
      "Error! Invalid decomp dimension.\n"
      "   - file name  : " + file.name + "\n"
      "   - file dims  : " + file_dims_str + "\n"
      "   - decomp dim : " + dim_name + "\n");
  EKAT_REQUIRE_MSG (file.decomps.find(dim_name)==file.decomps.end(),
      "Error! Decomposition was already added.\n"
      "   - file name  : " + file.name + "\n"
      "   - decomp dim : " + dim_name + "\n");

  auto dim = dim_it->second;
  int gcount;
  comm.all_reduce(&count,&gcount,1,MPI_SUM);

  EKAT_REQUIRE_MSG (count>=0 && gcount==dim->len,
      "Error! Invalid local start/count for decomposition.\n"
      "   - file name    : " + file.name + "\n"
      "   - decomp dim   : " + dim_name + "\n"
      "   - rank         : " + std::to_string(rank) + "\n"
      "   - local count  : " + std::to_string(count) + "\n"
      "   - global count : " + std::to_string(gcount) + "\n"
      "   - dim length   : " + std::to_string(dim->len) + "\n");

  std::vector<int> end (comm.size());
  end[rank] = count;
  comm.scan(end.data(),1,MPI_SUM);

  auto decomp = std::make_shared<NCDecomp>();
  decomp->dim = file.dims.at(dim_name);
  decomp->start = end[rank] - count;
  decomp->count = count;

  file.decomps[dim_name] = decomp;
}

bool has_var (const NCFile& file,const std::string& vname)
{
  return file.vars.find(vname)!=file.vars.end();
}

// template<>
// void write_var<int> (const NCFile& file,const std::string& vname,
//                      const int* const  data)
// {
//   const auto& var = file.vars.at(vname);
//   EKAT_REQUIRE_MSG (var.dtype=="int",
//       "Error! Input data type differs from declared var data.\n"
//       "   - file name   : " + file.name + "\n"
//       "   - var name    : " + vname + "\n"
//       "   - var dtype   : " + var.dtype + "\n"
//       "   - input dtype : int\n");

// }
// template<typename T>
// void read_var (const NCFile& file,
//                const std::string& vname,
//                T* const data);

} // namespace cldera
