#ifndef CLCDERA_PNETCDF_HPP
#define CLCDERA_PNETCDF_HPP

#include <ekat/mpi/ekat_comm.hpp>

#include <string>
#include <vector>
#include <memory>
#include <map>

namespace cldera {

// =============== TYPES =============== //

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

// A handy short name for string->T maps.
template<typename T>
using strmap_t = std::map<std::string,T>;


// A small struct to hold Netcdf dimensions info
struct NCDim {
  int           ncid;
  std::string   name;
  int           len;

  // Information for MPI-based decomposition (if any)
  bool          decomp = false;
  std::vector<int>  entries; // On-proc entries along partitioned dim
};

// A small struct to hold Netcdf variables info
struct NCVar {
  using dim_ptr_t = std::shared_ptr<const NCDim>;

  int                     ncid;
  std::string             name;
  std::string             dtype;
  std::vector<dim_ptr_t>  dims;
  void*                   data;
  int                     nrecords; // For time-dep vars only

  // If dims are not decomposed, each read/write handles the whole array,
  // so chunk_len is the prod of all (non-time) dim lengths.
  // If the 1st non-time dim is decomposed, each read/write handles
  // a single entry along the partitioned dim, so chunk_len is the prod of
  // all non-time (and non-decomposed) dim lengths.
  // Note: only the 1st non-time dimension can be decomposed.
  int                     chunk_len;
};

// A small struct to hold Netcdf file info
struct NCFile {
  using dim_ptr_t    = std::shared_ptr<NCDim>;
  using var_ptr_t    = std::shared_ptr<NCVar>;

  int         ncid;
  std::string name;
  IOMode      mode;

  strmap_t<dim_ptr_t>  dims;
  strmap_t<var_ptr_t>  vars;

  std::map<int,std::string> dim_id2name;

  ekat::Comm    comm;
};

// ================= FUNCTIONS ================ //

// --- File operations
std::shared_ptr<NCFile>
open_file (const std::string& fname, const ekat::Comm& comm, const IOMode mode);
void close_file (NCFile& file);

void add_dim (NCFile& file, const std::string& dname, const int len);

void add_var (      NCFile& file,
              const std::string& vname,
              const std::string& dtype,
                    std::vector<std::string> dims,
              const bool time_dep);

void add_decomp (      NCFile& file,
                 const std::string& dim_name,
                 const std::vector<int>& entries);

void enddef (const NCFile& file);
void redef  (const NCFile& file);

// --- Read/write operations
template<typename T>
void write_var (const NCFile& file,
                const std::string& vname,
                const T* const  data);

template<typename T>
void read_var (const NCFile& file,
               const std::string& vname,
                     T* const data,
               const int record = -1);

} // namespace cldera

#endif // CLCDERA_PNETCDF_HPP
