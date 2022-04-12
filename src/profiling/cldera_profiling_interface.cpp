#include "cldera_profiling_interface.hpp"

#include "cldera_profiling_session.hpp"
#include "cldera_profiling_archive.hpp"
#include "stats/cldera_compute_stat.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_parse_yaml_file.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

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

  s.clean_up();
}

void cldera_add_field_c (const char* name,
                         const int  rank,
                         const int* dims)
{
  cldera_add_partitioned_field_c(name,rank,dims,1,0);
}

void cldera_add_partitioned_field_c (
    const char* name,
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
    const char* name,
    const int   part,
    const int   part_size,
    const cldera::Real* data)
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  archive.get_field(name).set_part_data (part,part_size,data);
}

void cldera_set_field_c (
    const char* name,
    const cldera::Real* data)
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  archive.get_field(name).set_data (data);
}

void cldera_commit_field_c (const char* name)
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  archive.get_field(name).commit();
}

void cldera_commit_all_fields_c ()
{
  auto& s = cldera::ProfilingSession::instance(true);
  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  archive.commit_all_fields();
}

void cldera_init_requests_c ()
{
  auto& s = cldera::ProfilingSession::instance(true);

  using requests_t = std::map<std::string,std::vector<cldera::StatType>>;
  auto& requests = s.create<requests_t>("requests");

  //TODO: make the filename configurable
  ekat::ParameterList params = ekat::parse_yaml_file("./cldera_profiling_requests.yaml");

  for (int i=0; i<params.get<int>("Number Of Requests"); ++i) {
    auto key = ekat::strint("Request",i);
    const auto& req_pl = params.sublist(key);

    const auto& fname = req_pl.get<std::string>("Field Name");
    auto& req_stats = requests[fname];
    for (auto stat_str : req_pl.get<std::vector<std::string>>("Compute Stats")) {
      req_stats.push_back(cldera::str2stat (stat_str));
    }
  }
}

void cldera_compute_stats_c (const cldera::Real time)
{
  auto& s = cldera::ProfilingSession::instance(true);

  using requests_t = std::map<std::string,std::vector<cldera::StatType>>;
  auto& requests = s.get<requests_t>("requests");

  auto& archive = s.get<cldera::ProfilingArchive>("archive");

  for (const auto& it : requests) {
    const auto& fname = it.first;
    const auto& stats = it.second;
    const auto& f = archive.get_field(fname);
    for (auto stat : stats) {
      auto& history = archive.get_stat_history(fname,stat);
      cldera::compute_stat(time,f,stat,history);
    }
  }
}

} // extern "C"
