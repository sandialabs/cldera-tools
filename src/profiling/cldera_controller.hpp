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

class Controller
{

  using vos_t = std::vector<std::string>;
  using vod_t = std::vector<double>;

public:

  Controller(const ekat::ParameterList& params);

  bool is_ctrl_time(const int ymd, const int tod);

  // recording controller history
  // void hist_init();
  std::vector<double> get_ctrl_hist_at(const int ymd);

  // calculating controller
  void operator()(const int ymd);

private:

  // controller definition
  TimeStamp m_start_time, m_end_time;
  int m_ctrl_interval_months;
  vos_t m_ctrl_names;
  vos_t m_ctrl_types;
  std::vector<vod_t> m_ctrl_gains;
  vos_t m_ctrl_fields;
  vos_t m_ctrl_stats;
  vod_t m_ctrl_refvals;
  vos_t m_ctrl_projs;
  vos_t m_ff_types;
  vod_t m_ff_vals;
  std::vector<vos_t> m_cnstr_defs;

  // injection definitions
  std::vector<std::string> m_inj_names;
  std::vector<double> m_inj_lats;
  std::vector<double> m_inj_lons;
  std::vector<double> m_inj_alts;
  std::vector<vod_t> m_inj_transfer_funcs;

  // controller history
  std::map<int, std::vector<double>> m_ctrl_hist;


};


}

#endif