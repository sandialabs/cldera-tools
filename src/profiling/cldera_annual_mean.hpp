#ifndef CLDERA_ANNUAL_MEAN_HPP
#define CLDERA_ANNUAL_MEAN_HPP

#include "cldera_time_stamp.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>
#include <string>
#include <vector>

namespace cldera {

/*
 * AnnualMeanAccumulator (BMW, G6-1.5K-HiLLA annual PI controller)
 *
 * Accumulates calendar-year means of scalar stats (e.g. TREFHT_glb, the global-mean
 * 2 m temperature) from EVERY cldera_compute_stats call, so that an annual feedback
 * controller can use the mean of the previous year.
 *
 * Why a restart file is needed:
 *   cldera_compute_stats is called with the time stamp at the START of each step. In a
 *   run that stops at Jan 1 00:00 (annual restarts), the last call is stamped Dec 31
 *   23:30 (dt=1800 s), so the year is complete only at the end of the run, and the
 *   annual decision is taken at the first call of the NEXT run. The accumulator state
 *   is therefore written at the end of every run (cldera_clean_up), to
 *       <prefix>.<stop time YYYY-MM-DD-SSSSS>.txt
 *   and read at the start of the next run, from <prefix>.<run start>.txt.
 *   Naming the file by the restart time (like E3SM restart files) means that re-running
 *   from an older restart picks up the matching accumulator file.
 *
 *   - If the run starts at the case start (startup/branch/hybrid), no file is read.
 *   - Otherwise the file for the run start time must exist (unless 'required' is false,
 *     in which case the accumulator starts empty, with a warning).
 *   - The file is checked for continuity: its last sample must be at most 2 steps
 *     before the first sample of this run (the run_t0 call itself is skipped by cldera).
 *
 * Only the root rank does any work: the stats are global (all-reduced) and the Python
 * controller only runs on root.
 *
 * Calendar: EAM uses a 365-day (noleap) calendar; all time differences use days_365.
 */
class AnnualMeanAccumulator {
public:
  struct StatAccum {
    std::string name;
    // Current (open) year
    int    year     = -1;
    long   nsamples = 0;
    double sum      = 0;
    // Last completed (closed) year
    int    closed_year     = -1;
    long   closed_nsamples = 0;
    double closed_mean     = std::numeric_limits<double>::quiet_NaN();
    double closed_coverage = std::numeric_limits<double>::quiet_NaN();
  };

  AnnualMeanAccumulator (const ekat::Comm& comm,
                         const std::vector<std::string>& stat_names,
                         const std::string& file_prefix,
                         const bool restart_required,
                         const TimeStamp& case_t0,
                         const TimeStamp& run_t0);

  // Add one sample of every stat (same order as stat_names). Closes the open year
  // if t is in a later year.
  void add (const TimeStamp& t, const std::vector<double>& vals);

  // Write the restart file for the run that starts at 'stop'.
  void write_restart (const TimeStamp& stop) const;

  const std::vector<std::string>& stat_names () const { return m_names; }
  const std::vector<StatAccum>& stats () const { return m_stats; }
  double dt_seconds () const { return m_dt; }
  bool   has_samples () const { return m_has_last; }
  const TimeStamp& last_time () const { return m_last; }
  const std::string& loaded_file () const { return m_loaded_file; }
  std::string restart_file_name (const TimeStamp& t) const;

private:
  void read_restart (const std::string& fname);

  ekat::Comm               m_comm;
  std::vector<std::string> m_names;
  std::vector<StatAccum>   m_stats;
  std::string              m_prefix;
  TimeStamp                m_case_t0, m_run_t0;

  double    m_dt = std::numeric_limits<double>::quiet_NaN();  // [s], smallest sample spacing seen
  TimeStamp m_last;            // time of the last sample
  bool      m_has_last = false;
  bool      m_first_add_of_run = true;
  std::string m_loaded_file;
};

} // namespace cldera

#endif // CLDERA_ANNUAL_MEAN_HPP
