#include "cldera_profiling_interface.hpp"

#include "cldera_profiling_session.hpp"
#include "cldera_profiling_archive.hpp"
#include "stats/cldera_compute_stat.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_parse_yaml_file.hpp>
#include <ekat/util/ekat_string_utils.hpp>

extern "C" {

void cldera_profiling_init (const MPI_Fint& fcomm)
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

void cldera_profiling_clean_up ()
{
  auto& s = cldera::ProfilingSession::instance();
  if (s.get_comm().am_i_root()) {
    printf(" -> Shutting down cldera profiling session...\n");
  }

  s.clean_up();
}

void cldera_add_field (const char* name,
                       const cldera::Real*& data,
                       const int& rank,
                       const int*& dims)
{
  auto& s = cldera::ProfilingSession::instance(true);

  // Copy input raw pointer to vector
  std::vector<int> d(rank);
  for (int i=0; i<rank; ++i) {
    d[i] = dims[i];
  }

  // Set data in the archive structure
  auto& archive = s.create_or_get<cldera::ProfilingArchive>("archive");
  archive.add_field(name,data,d);
}

void cldera_init_requests ()
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

void cldera_compute_stats (const cldera::Real& time)
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
