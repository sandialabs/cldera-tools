#include "cldera_annual_mean.hpp"

#include <ekat/ekat_assert.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace cldera {

namespace {
constexpr double kSecPerYear365 = 365.0*86400.0;

double secs_between (const TimeStamp& a, const TimeStamp& b) {
  // b - a [s], noleap calendar
  return (days_365(b.ymd(),b.tod()) - days_365(a.ymd(),a.tod()))*86400.0;
}
} // anonymous namespace

AnnualMeanAccumulator::
AnnualMeanAccumulator (const ekat::Comm& comm,
                       const std::vector<std::string>& stat_names,
                       const std::string& file_prefix,
                       const bool restart_required,
                       const TimeStamp& case_t0,
                       const TimeStamp& run_t0)
 : m_comm (comm)
 , m_names (stat_names)
 , m_prefix (file_prefix)
 , m_case_t0 (case_t0)
 , m_run_t0 (run_t0)
{
  for (const auto& n : m_names) {
    StatAccum s;
    s.name = n;
    m_stats.push_back(s);
  }

  // Decide collectively whether a restart file is needed and present, so that ALL ranks
  // error out together (only root reads the file).
  const bool at_case_start = (m_run_t0==m_case_t0);
  const auto fname = restart_file_name(m_run_t0);
  int found = 0;
  if (m_comm.am_i_root() and not at_case_start) {
    found = std::ifstream(fname).good() ? 1 : 0;
  }
  m_comm.broadcast(&found,1,m_comm.root_rank());
  EKAT_REQUIRE_MSG (at_case_start or found==1 or not restart_required,
      "[AnnualMean] Error! This run does not start at the case start, but the annual-mean\n"
      "  restart file for the run start time was not found.\n"
      " - expected file: " + fname + "\n"
      " - case start   : " + m_case_t0.to_string() + "\n"
      " - run start    : " + m_run_t0.to_string() + "\n"
      "  The file is written at the end of every run and at every Jan 1. If you restarted\n"
      "  from a restart time that is neither, the annual mean of the current year is lost.\n"
      "  Set 'Annual Mean Restart Required: false' to start with no history (the current\n"
      "  year's mean will then be incomplete and the controller will refuse to use it).\n");

  if (not m_comm.am_i_root()) return;

  if (at_case_start) {
    printf(" [CLDERA] AnnualMean: run starts at the case start (%s); starting with no history.\n",
           m_run_t0.to_string().c_str());
    return;
  }
  if (found==1) {
    read_restart(fname);
  } else {
    printf(" [CLDERA] AnnualMean WARNING: restart file '%s' not found; starting with no history.\n",
           fname.c_str());
  }
}

std::string AnnualMeanAccumulator::
restart_file_name (const TimeStamp& t) const
{
  return m_prefix + "." + t.to_string() + ".txt";
}

void AnnualMeanAccumulator::
add (const TimeStamp& t, const std::vector<double>& vals)
{
  if (not m_comm.am_i_root()) return;
  EKAT_REQUIRE_MSG (vals.size()==m_stats.size(),
      "[AnnualMean] Error! Wrong number of values.\n");

  if (m_has_last) {
    const double ds = secs_between(m_last,t);
    EKAT_REQUIRE_MSG (ds>0,
        "[AnnualMean] Error! Samples must be strictly increasing in time.\n"
        " - last sample: " + m_last.to_string() + "\n"
        " - this sample: " + t.to_string() + "\n");
    if (m_first_add_of_run) {
      // Continuity across a restart: the run_t0 call is skipped by cldera, so the first
      // sample of a run is normally 2 steps after the last sample of the previous run.
      if (std::isfinite(m_dt)) {
        EKAT_REQUIRE_MSG (ds <= 2.0*m_dt + 1.0,
            "[AnnualMean] Error! Gap between the restart file and this run.\n"
            " - last sample in restart file: " + m_last.to_string() + "\n"
            " - first sample of this run   : " + t.to_string() + "\n"
            " - time step                  : " + std::to_string(m_dt) + " s\n"
            "  Was the restart file written by a different run?\n");
      }
    } else {
      m_dt = std::isfinite(m_dt) ? std::min(m_dt,ds) : ds;
    }
  }
  m_first_add_of_run = false;

  const int y = t.year();
  bool closed_now = false;
  for (size_t i=0; i<m_stats.size(); ++i) {
    auto& s = m_stats[i];
    if (s.year>=0 and y>s.year) {
      // Close the open year
      s.closed_year     = s.year;
      s.closed_nsamples = s.nsamples;
      s.closed_mean     = s.nsamples>0 ? s.sum/s.nsamples : std::numeric_limits<double>::quiet_NaN();
      s.closed_coverage = std::isfinite(m_dt) ? s.nsamples*m_dt/kSecPerYear365
                                              : std::numeric_limits<double>::quiet_NaN();
      printf(" [CLDERA] AnnualMean: closed year %04d for %s: mean=%.6f, nsamples=%ld, coverage=%.5f\n",
             s.closed_year, s.name.c_str(), s.closed_mean, s.closed_nsamples, s.closed_coverage);
      s.year = -1;
      closed_now = true;
    }
  }
  // If a year was just closed INSIDE a run, also write a restart file for Jan 1 00:00 of the
  // new year, with the state BEFORE this sample. This supports restarting from any Jan 1
  // restart, not only from the end of a run. Skip it if that file is the one we just read.
  if (closed_now) {
    const TimeStamp jan1(y*10000+101,0);
    if (restart_file_name(jan1)!=m_loaded_file) {
      write_restart(jan1);
    }
  }

  for (size_t i=0; i<m_stats.size(); ++i) {
    auto& s = m_stats[i];
    if (s.year<0) {
      s.year = y;
      s.nsamples = 0;
      s.sum = 0;
    }
    EKAT_REQUIRE_MSG (std::isfinite(vals[i]),
        "[AnnualMean] Error! Non-finite value for stat '" + s.name + "' at " + t.to_string() + "\n");
    s.sum += vals[i];
    ++s.nsamples;
  }
  m_last = t;
  m_has_last = true;
}

void AnnualMeanAccumulator::
write_restart (const TimeStamp& stop) const
{
  if (not m_comm.am_i_root()) return;
  const auto fname = restart_file_name(stop);
  const auto tmp = fname + ".tmp";
  {
    std::ofstream out(tmp);
    EKAT_REQUIRE_MSG (out.good(), "[AnnualMean] Error! Could not open '" + tmp + "' for writing.\n");
    out << std::setprecision(17);
    out << "# cldera-tools AnnualMeanAccumulator restart file (BMW). Do not edit.\n";
    out << "version 1\n";
    out << "restart_time " << stop.ymd() << " " << stop.tod() << "\n";
    out << "case_t0 " << m_case_t0.ymd() << " " << m_case_t0.tod() << "\n";
    out << "dt " << (std::isfinite(m_dt) ? m_dt : -1.0) << "\n";
    out << "last " << (m_has_last ? 1 : 0) << " " << m_last.ymd() << " " << m_last.tod() << "\n";
    out << "nstats " << m_stats.size() << "\n";
    for (const auto& s : m_stats) {
      out << "stat " << s.name << " "
          << s.year << " " << s.nsamples << " " << s.sum << " "
          << s.closed_year << " " << s.closed_nsamples << " "
          << (std::isfinite(s.closed_mean) ? s.closed_mean : -1.0e300) << " "
          << (std::isfinite(s.closed_coverage) ? s.closed_coverage : -1.0) << "\n";
    }
  }
  std::rename(tmp.c_str(),fname.c_str());
  printf(" [CLDERA] AnnualMean: wrote restart file '%s'\n", fname.c_str());
}

void AnnualMeanAccumulator::
read_restart (const std::string& fname)
{
  std::ifstream in(fname);
  std::string line, key;
  int version = -1;
  std::vector<StatAccum> loaded;
  while (std::getline(in,line)) {
    if (line.empty() or line[0]=='#') continue;
    std::istringstream ss(line);
    ss >> key;
    if (key=="version") {
      ss >> version;
    } else if (key=="restart_time") {
      int ymd, tod;
      ss >> ymd >> tod;
      EKAT_REQUIRE_MSG (TimeStamp(ymd,tod)==m_run_t0,
          "[AnnualMean] Error! Restart file time does not match the run start.\n"
          " - file: " + fname + "\n");
    } else if (key=="case_t0") {
      int ymd, tod;
      ss >> ymd >> tod;
      EKAT_REQUIRE_MSG (TimeStamp(ymd,tod)==m_case_t0,
          "[AnnualMean] Error! Restart file was written by a case with a different start time.\n"
          " - file: " + fname + "\n");
    } else if (key=="dt") {
      double dt;
      ss >> dt;
      m_dt = dt>0 ? dt : std::numeric_limits<double>::quiet_NaN();
    } else if (key=="last") {
      int has, ymd, tod;
      ss >> has >> ymd >> tod;
      m_has_last = has==1;
      if (m_has_last) m_last = TimeStamp(ymd,tod);
    } else if (key=="stat") {
      StatAccum s;
      double cm, cc;
      ss >> s.name >> s.year >> s.nsamples >> s.sum >> s.closed_year >> s.closed_nsamples >> cm >> cc;
      s.closed_mean = cm>-1.0e299 ? cm : std::numeric_limits<double>::quiet_NaN();
      s.closed_coverage = cc>=0 ? cc : std::numeric_limits<double>::quiet_NaN();
      loaded.push_back(s);
    }
  }
  EKAT_REQUIRE_MSG (version==1, "[AnnualMean] Error! Unsupported restart file version in '" + fname + "'.\n");
  for (auto& s : m_stats) {
    bool found = false;
    for (const auto& l : loaded) {
      if (l.name==s.name) { s = l; found = true; }
    }
    EKAT_REQUIRE_MSG (found,
        "[AnnualMean] Error! Stat '" + s.name + "' not found in restart file '" + fname + "'.\n");
  }
  m_loaded_file = fname;
  printf(" [CLDERA] AnnualMean: read restart file '%s' (last sample %s, dt=%g s)\n",
         fname.c_str(), m_last.to_string().c_str(), m_dt);
  for (const auto& s : m_stats) {
    printf(" [CLDERA]   %s: open year %04d (%ld samples), last closed year %04d (mean=%.6f)\n",
           s.name.c_str(), s.year, s.nsamples, s.closed_year, s.closed_mean);
  }
}

} // namespace cldera
