#ifndef CLDERA_CONTROLLER_HPP
#define CLDERA_CONTROLLER_HPP

#include <string>
#include <vector>
#include <limits>
#include <map>
#include <mpi.h>

#include "ekat/io/ekat_yaml.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include "cldera_time_stamp.hpp"

namespace cldera {

class ProfilingContext;

struct IntervalDef {
  std::string name;
  int start_mmdd;
  int end_mmdd;
  int span_days;
};

struct ControllerDef {
  std::string name;      // controller name (e.g., T0)
  std::string interval;  // interval name (e.g., DJF)
  std::string type;

  // constant controller
  double const_val = std::numeric_limits<double>::quiet_NaN();

  // non-constant controller
  std::vector<double> gains;
  std::string field;
  std::string stat;
  double refval = std::numeric_limits<double>::quiet_NaN();
  std::string proj;
  std::string ff_type;
  double ff_val = std::numeric_limits<double>::quiet_NaN();
  std::vector<std::string> cnstr_defs;

  // BMW (G6-1.5K-HiLLA): annual PI feedback controller options
  //  - has_refval: 'Reference' was given in the yaml. If not, the Python controller
  //    computes it from the annual mean of the year BEFORE the controller starts,
  //    plus ref_offset ('Reference Offset', default -1.0, i.e. 1 K below that year).
  //  - init_annual_mean: 'Initial Annual Mean'. Annual mean of the year before the
  //    controller starts, if it is not simulated (e.g. a branch run from a control run).
  //  - ff_start_year: 'Feedforward: Start Year'. Year with FF = 1*ff_val (linear FF).
  //    Default (-1): the controller start year.
  //  - fb_start_year: 'Feedback Start Year'. First year whose decision uses feedback
  //    (P/I terms, from the annual mean of the year before). In earlier years of the
  //    controller window the output is the feedforward only, and no error is recorded.
  //    Default (-1): the controller start year (feedback from the first decision).
  //  - min_out/max_out: 'Minimum Output'/'Maximum Output'. The controller output
  //    (requested global-mean AOD) is clipped to [min_out, max_out]. Default [0, none].
  //  - min_coverage: 'Minimum Annual Coverage'. Minimum fraction of a year that must be
  //    sampled for its annual mean to be used. Default 0.99.
  bool   has_refval = false;
  double ref_offset = -1.0;
  double init_annual_mean = std::numeric_limits<double>::quiet_NaN();
  int    ff_start_year = -1;
  int    fb_start_year = -1;
  double min_out = 0.0;
  double max_out = std::numeric_limits<double>::quiet_NaN();
  double min_coverage = 0.99;
};

class Controller
{

  using vos_t = std::vector<std::string>;
  using vod_t = std::vector<double>;
  using voi_t = std::vector<int>;

public:

  // case_t0 (the case start time) is needed only if the controller yaml uses
  // 'Start Years After Case Start' instead of 'Start Day'.
  Controller(const ekat::ParameterList& params);
  Controller(const ekat::ParameterList& params, const TimeStamp& case_t0);

  bool is_ctrl_time(const int ymd, const int tod);
  bool is_interval_start(const int ymd, const int tod) const;

  // recording controller history
  // void hist_init();
  vod_t get_ctrl_hist_at(const int ymd);
  bool has_ctrl_hist_at(const int ymd) const { return m_ctrl_hist.find(ymd)!=m_ctrl_hist.end(); }
  // Overwrite history entry (used to store the broadcast result on ALL ranks)
  void set_ctrl_hist_at(const int ymd, const vod_t& v) { m_ctrl_hist[ymd] = v; }

  // calculating controller
  vod_t operator()(const int ymd, const int tod);

  // control intervals
  const std::vector<IntervalDef>& get_intervals() const { return m_intervals; }

  // getting controller definitions
  const std::vector<ControllerDef>& get_controllers() const { return m_controllers; }

  // Unique names of the stats used by the feedback (non-constant) controllers, in order of
  // appearance. cldera-tools accumulates their annual means.
  vos_t get_feedback_stats() const;

  // getting injection definitions
  const vos_t get_inj_names() { return m_inj_names; }
  const vod_t get_inj_lats() { return m_inj_lats; }
  const vod_t get_inj_lons() { return m_inj_lons; }
  const vod_t get_inj_alts() { return m_inj_alts; }
  const std::vector<vod_t> get_inj_transfer_funcs() {return m_inj_transfer_funcs; }

  // Controller active window (resolved), and case start (if known)
  int get_start_ymd() const { return m_start_time.ymd(); }
  int get_start_tod() const { return m_start_time.tod(); }
  int get_end_ymd()   const { return m_end_time.ymd(); }
  int get_end_tod()   const { return m_end_time.tod(); }
  int get_case_start_ymd() const { return m_case_t0.ymd(); }
  int get_case_start_tod() const { return m_case_t0.tod(); }
  bool has_case_start() const { return m_has_case_t0; }
  // How the start was specified: "Start Day" or "Start Years After Case Start"
  const std::string& get_start_mode() const { return m_start_mode; }
  int get_start_years_after_case_start() const { return m_start_years_after; }

  // Name of the controller yaml file
  const std::string& get_filename() const { return m_filename; }

  void set_context(ProfilingContext* ctx) { m_context = ctx; }

  // Info about the current call, passed to Python (first call of this run, run start)
  void set_call_info(const bool first_call_of_run, const TimeStamp& run_t0) {
    m_first_call_of_run = first_call_of_run;
    m_run_t0 = run_t0;
  }
  bool is_first_call_of_run() const { return m_first_call_of_run; }
  int get_run_start_ymd() const { return m_run_t0.ymd(); }
  int get_run_start_tod() const { return m_run_t0.tod(); }

private:

  void init (const ekat::ParameterList& params);

  // controller definition
  TimeStamp m_start_time, m_end_time;
  TimeStamp m_case_t0;
  bool      m_has_case_t0 = false;
  std::string m_start_mode;
  int       m_start_years_after = -1;
  std::string m_filename;
  std::vector<IntervalDef>   m_intervals;
  std::vector<ControllerDef> m_controllers;

  // injection definitions
  vos_t m_inj_names;
  vod_t m_inj_lats;
  vod_t m_inj_lons;
  vod_t m_inj_alts;
  std::vector<vod_t> m_inj_transfer_funcs;

  // controller history
  std::map<int, std::vector<double>> m_ctrl_hist;

  ProfilingContext* m_context = nullptr;

  bool      m_first_call_of_run = false;
  TimeStamp m_run_t0;

};


}

#endif
