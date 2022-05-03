#include "cldera_profiling_interface.hpp"

#include "cldera_profiling_session.hpp"
#include "cldera_profiling_archive.hpp"
#include "stats/cldera_compute_stat.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <cstring>
#include <fstream>

extern "C" {

void cldera_init_c (const MPI_Fint fcomm)
{
  // Convert F90 comm to C comm, create ekat::Comm, and init session
  MPI_Comm mpiComm = MPI_Comm_f2c(fcomm);
  EKAT_REQUIRE_MSG (mpiComm!=nullptr, "Error! Input fortran comm seems invalid.\n");

  ekat::Comm comm(mpiComm);

  if (comm.am_i_root()) {
    printf(" -> Initializing cldera profiling session...\n");
  }

  auto& s = cldera::ProfilingSession::instance();
  s.init(comm);

  using requests_t = std::map<std::string,std::vector<cldera::StatType>>;
  using vos_t = std::vector<std::string>;

  //TODO: make the filename configurable
  std::string filename = "./cldera_profiling_config.yaml";
  if (std::ifstream(filename).good()) {
    auto& params = s.create<ekat::ParameterList>("params",ekat::parse_yaml_file(filename));

    auto& requests = s.create<requests_t>("requests");
    const auto& fnames = params.get<vos_t>("Fields To Track");
    for (const auto& fname : fnames) {
      const auto& req_pl = params.sublist(fname);
      auto& req_stats = requests[fname];
      for (auto stat_str : req_pl.get<vos_t>("Compute Stats")) {
        req_stats.push_back(cldera::str2stat (stat_str));
      }
    }
  } else {
    if (comm.am_i_root()) {
      printf(" -> WARNING: no 'cdlera_profiling_config.yaml' file found.\n"
             "    CLDERA profiling tools will do nothing.\n");
    }

    // Create entries, since they will be queried later.
    // Since they're empty, all ops will be trivial/no-ops.
    auto& params = s.create<ekat::ParameterList>("params");
    params.set<std::vector<std::string>>("Fields To Track",{});
    s.create<requests_t>("requests");
  }

  if (comm.am_i_root()) {
    printf(" -> Initializing cldera profiling session...done!\n");
  }
}

void cldera_clean_up_c ()
{
  auto& s = cldera::ProfilingSession::instance();
  if (s.get_comm().am_i_root()) {
    printf(" -> Shutting down cldera profiling session...\n");
  }

  auto& archive = s.create_or_get<cldera::ProfilingArchive>("archive");
  const auto& p = s.get<ekat::ParameterList>("params");
  if (p.isParameter("Stats Output File")) {
    const auto& ofile = p.get<std::string>("Stats Output File");
    archive.dump_stats_to_yaml (ofile);
  }

  s.clean_up();
  if (s.get_comm().am_i_root()) {
    printf(" -> Shutting down cldera profiling session...done!\n");
  }
}

void cldera_add_field_c (const char*& name,
                         const int  rank,
                         const int* dims)
{
  cldera_add_partitioned_field_c(name,rank,dims,1,0);
}

void cldera_add_partitioned_field_c (
    const char*& name,
    const int  rank,
    const int* dims,
    const int  num_parts,
    const int  part_dim)
{
  EKAT_REQUIRE_MSG (rank>=0 && rank<=4,
      "Error! Unsupported field rank (" + std::to_string(rank) + "\n");
  EKAT_REQUIRE_MSG (num_parts>=1,
      "Error! Invalid number of partitions (" + std::to_string(num_parts) + "\n");

  auto& s = cldera::ProfilingSession::instance(true);

  // Copy input raw pointer to vector
  std::vector<int> d(rank);
  for (int i=0; i<rank; ++i) {
    EKAT_REQUIRE_MSG (dims[i]>=0,
        "Error! Invalid field extent.\n"
        "   - Field name:" +  std::string(name) + "\n";
        "   - Dimension: " +  std::to_string(i) + "\n"
        "   - Extent:    " +  std::to_string(dims[i]) + "\n");

    d[i] = dims[i];
  }

  // Set data in the archive structure
  auto& archive = s.create_or_get<cldera::ProfilingArchive>("archive");
  archive.add_field(cldera::Field(name,d,num_parts,part_dim));
}

void cldera_set_field_partition_c (
    const char*& name,
    const int   part,
    const int   part_size,
    const cldera::Real*& data)
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  archive.get_field(name).set_part_data (part,part_size,data);
}

void cldera_set_field_c (
    const char*& name,
    const cldera::Real*& data)
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  archive.get_field(name).set_data (data);
}

void cldera_commit_field_c (const char*& name)
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  archive.get_field(name).commit();
}

void cldera_commit_all_fields_c ()
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.create_or_get<cldera::ProfilingArchive>("archive");

  archive.commit_all_fields();
}

int cldera_get_num_fields_c ()
{
  const auto& s = cldera::ProfilingSession::instance(true);
  const auto& p = s.get<ekat::ParameterList>("params");
  return p.get<std::vector<std::string>>("Fields To Track").size();
}

void cldera_get_field_name_c (const int i, char*& name)
{
  const auto& s = cldera::ProfilingSession::instance(true);
  const auto& p = s.get<ekat::ParameterList>("params");
  const auto& names = p.get<std::vector<std::string>>("Fields To Track");

  EKAT_REQUIRE_MSG (i>=0 && i<cldera_get_num_fields_c(),
      "Error! Invalid field index.\n"
      "  - ifield: " + std::to_string(i) + "\n"
      "  - num fields: " + std::to_string(cldera_get_num_fields_c()) + "\n");
  strcpy(name,names[i].data());
}

void cldera_compute_stats_c (const int ymd, const int tod)
{
  auto& s = cldera::ProfilingSession::instance(true);
  const auto& comm = s.get_comm();

  using requests_t = std::map<std::string,std::vector<cldera::StatType>>;
  auto& requests = s.get<requests_t>("requests");

  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  for (const auto& it : requests) {
    const auto& fname = it.first;
    const auto& stats = it.second;
    const auto& f = archive.get_field(fname);
    for (auto stat : stats) {
      auto& history = archive.get_stat_history(fname,stat);
      cldera::compute_stat({ymd,tod},f,stat,history,comm);
    }
  }
}

} // extern "C"
