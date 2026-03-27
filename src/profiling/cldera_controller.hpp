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
};

class Controller
{

  using vos_t = std::vector<std::string>;
  using vod_t = std::vector<double>;
  using voi_t = std::vector<int>;

public:

  Controller(const ekat::ParameterList& params);

  bool is_ctrl_time(const int ymd, const int tod);

  // recording controller history
  // void hist_init();
  vod_t get_ctrl_hist_at(const int ymd);

  // calculating controller
  vod_t operator()(const int ymd, const int tod);

  // control intervals
  const std::vector<IntervalDef>& get_intervals() const { return m_intervals; }

  // getting controller definitions
  const std::vector<ControllerDef>& get_controllers() const { return m_controllers; }

  // getting injection definitions
  const vos_t get_inj_names() { return m_inj_names; }
  const vod_t get_inj_lats() { return m_inj_lats; }
  const vod_t get_inj_lons() { return m_inj_lons; }
  const vod_t get_inj_alts() { return m_inj_alts; }
  const std::vector<vod_t> get_inj_transfer_funcs() {return m_inj_transfer_funcs; }

  void set_context(ProfilingContext* ctx) { m_context = ctx; }

private:

  // controller definition
  TimeStamp m_start_time, m_end_time;
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

};


}

#endif
