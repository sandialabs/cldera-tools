#ifndef CLCDERA_PNETCDF_HPP
#define CLCDERA_PNETCDF_HPP

#include <ekat/mpi/ekat_comm.hpp>

#include <string>
#include <vector>
#include <memory>
#include <map>

namespace cldera {

// =============== TYPES =============== //

// A handy short name for string->T maps.
template<typename T>
using strmap_t = std::map<std::string,T>;

// A small struct to hold Netcdf dimensions info
struct NCDim {
  int           ncid;
  std::string   name;
  int           len;
};

// A small struct to hold Netcdf variables info
struct NCVar {
  using dim_ptr_t = std::shared_ptr<const NCDim>;

  int                     ncid;
  std::string             name;
  std::string             dtype;
  std::vector<dim_ptr_t>  dims;
};

// A small struct to hold Netcdf decomposition info
struct NCDecomp {
  using dim_ptr_t = std::shared_ptr<const NCDim>;

  dim_ptr_t     dim;
  int           start;
  int           count;
};

// An enum is a bit more verbose than int/bool flags
enum class IOMode {
  Read,
  Write,
  Invalid
};

inline std::string e2str (const IOMode m) {
  std::string name;
  switch (m) {
    case IOMode::Read:    name = "read";    break;
    case IOMode::Write:   name = "write";   break;
    case IOMode::Invalid: name = "invalid"; break;
  }
  return name;
}

// A small struct to hold Netcdf file info
struct NCFile {
  using dim_ptr_t    = std::shared_ptr<NCDim>;
  using var_ptr_t    = std::shared_ptr<NCVar>;
  using decomp_ptr_t = std::shared_ptr<NCDecomp>;

  int         ncid;
  std::string name;
  IOMode      mode;

  strmap_t<dim_ptr_t>  dims;
  strmap_t<var_ptr_t>  vars;

  std::map<int,std::string> dim_id2name;

  strmap_t<decomp_ptr_t>  decomps;

  ekat::Comm    comm;
};

// ================= FUNCTIONS ================ //

// File operations
std::shared_ptr<NCFile>
open_file (const std::string& fname, const ekat::Comm& comm, const IOMode mode);
void close_file (NCFile& file);
void enddef (const NCFile& file);
void redef  (const NCFile& file);

// Dimension operations
void add_dim (NCFile& file, const std::string& dname, const int len);
bool has_dim (const NCFile& file, const std::string& dname);
int get_dim_len (const NCFile& file, const std::string& dname);

// Variable operations
void add_var (      NCFile& file,
              const std::string& vname,
              const std::string& dtype,
              const std::vector<std::string>& dims);

void add_decomp (      NCFile& file,
                 const std::string& dim_name,
                 const int count);

bool has_var (const NCFile& file,const std::string& vname);

template<typename T>
void write_var (const NCFile& file,const std::string& vname,
                const T* const  data);
template<typename T>
void read_var (const NCFile& file,
               const std::string& vname,
               T* const data);

} // namespace cldera

#endif // CLCDERA_PNETCDF_HPP
