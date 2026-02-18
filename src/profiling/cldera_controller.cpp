
#include <iostream>

#include "cldera_controller.hpp"

namespace cldera {

Controller::
Controller(const ekat::ParameterList& params)
{

  using vos_t = std::vector<std::string>;
  using vod_t = std::vector<double>;

	// load controller parameters
  const auto prof_fields = params.get<vos_t>("Fields To Track");
	const auto ctrl_params_filename = params.get<std::string>("Controller");
  const auto ctrl_params = ekat::parse_yaml_file(ctrl_params_filename);

	// get time stamps
	const int start_ymd = ctrl_params.get<int>("Start Day");
	const int start_tod = ctrl_params.get<int>("Start Time");
	const int end_ymd = ctrl_params.get<int>("End Day");
	const int end_tod = ctrl_params.get<int>("End Time");
	m_start_time = TimeStamp(start_ymd, start_tod);
	m_end_time = TimeStamp(end_ymd, end_tod);

	// get controller definitions
  m_ctrl_interval_months = ctrl_params.get<int>("Interval Months");
  EKAT_REQUIRE_MSG(
    m_ctrl_interval_months >= 1,
    "Interval Months must be at least one month"
  );

	m_ctrl_names = ctrl_params.get<vos_t>("Controllers");
  size_t nctrl = m_ctrl_names.size();
	for (std::string& ctrl_name : m_ctrl_names) {

    auto ctrl_def = ctrl_params.sublist(ctrl_name);

    // type
    // TODO: this needs to be updated if using MPC
    auto ctrl_type = ctrl_def.get<std::string>("Type");
    EKAT_REQUIRE_MSG(
      ctrl_type == "P" || ctrl_type == "PI" || ctrl_type == "PD" || ctrl_type == "PID",
      "Invalid controller type for definition \"" + ctrl_name + "\" must be one of: P, PI, PD, PID"
    );
    m_ctrl_types.push_back(ctrl_type);

    // gains
    // TODO: need to be updated if using MPC
    auto ctrl_gains = ctrl_def.get<vod_t>("Gains");
    EKAT_REQUIRE_MSG(
      ctrl_gains.size() == ctrl_type.size(),
      "Incorrect number of gains for definition \"" + ctrl_name
    );
    m_ctrl_gains.push_back(ctrl_gains);

    // statistic
    auto ctrl_field = ctrl_def.get<std::string>("Field");
    EKAT_REQUIRE_MSG(
      std::find(prof_fields.begin(), prof_fields.end(), ctrl_field) != prof_fields.end(),
      "Controller field \"" + ctrl_field + "\" not found in tracked fields"
    );
    auto field_def = params.sublist(ctrl_field);
    auto field_stats = field_def.get<vos_t>("Compute Stats");
    auto ctrl_stat = ctrl_def.get<std::string>("Stat");
    EKAT_REQUIRE_MSG(
      std::find(field_stats.begin(), field_stats.end(), ctrl_stat) != field_stats.end(),
      "Invalid controller statistic \"" << ctrl_stat << "\" for tracked field \"" << ctrl_field << "\"" << std::endl;
    );
    m_ctrl_fields.push_back(ctrl_field);
    m_ctrl_stats.push_back(ctrl_stat);

    // residual reference value
    auto ctrl_refval = ctrl_def.get<double>("Reference");
    m_ctrl_refvals.push_back(ctrl_refval);

    // TODO: some checking on projections
    auto ctrl_proj = ctrl_def.get<std::string>("Projection");
    m_ctrl_projs.push_back(ctrl_proj);

    // TODO: this needs to be generalized greatly
    if (ctrl_def.isSublist("Feedforward")) {
      auto ff_def = ctrl_def.sublist("Feedforward");
      auto ff_type = ff_def.get<std::string>("Type");
      auto ff_val = ff_def.get<double>("Value");
      m_ff_types.push_back(ff_type);
      m_ff_vals.push_back(ff_val);
    } else {
      m_ff_types.push_back("");
      m_ff_vals.push_back(
        std::numeric_limits<double>::quiet_NaN()
      );
    }

    // constraints
    if (ctrl_def.isParameter("Constraints")) {
      auto cnstr_names = ctrl_def.get<vos_t>("Constraints");
      vos_t cnstr_defs;
      for (std::string& cnstr_name : cnstr_names) {
        auto cnstr_def = ctrl_def.get<std::string>("cnstr_" + cnstr_name);
        cnstr_defs.push_back(cnstr_def);
      }
      m_cnstr_defs.push_back(cnstr_defs);
    }

  }

  // injection definitions
  m_inj_names = ctrl_params.get<vos_t>("Injections");
  auto tf_defs = ctrl_params.sublist("Transfer Function");
	for (std::string& inj_name : m_inj_names) {

    auto inj_def = ctrl_params.sublist(inj_name);

    auto lat = inj_def.get<double>("lat");
    EKAT_REQUIRE_MSG(
      lat >= -90.0 && lat <= 90,
      "Invalid injection latitude for injection label \"" + inj_name << "\": " + std::to_string(lat)
    );
    m_inj_lats.push_back(lat);

    auto lon = inj_def.get<double>("lon");
    EKAT_REQUIRE_MSG(
      lon >= 0.0 && lon <= 360.0,
      "Invalid injection longitude for injection label \"" + inj_name << "\": " + std::to_string(lon)
    );
    m_inj_lons.push_back(lon);

    auto alt = inj_def.get<double>("alt");
    EKAT_REQUIRE_MSG(
      alt >= 0.0,
      "Invalid injection altitude for injection label \"" + inj_name << "\": " + std::to_string(alt)
    );
    m_inj_alts.push_back(alt);

    // transfer function
    auto tf = tf_defs.get<vod_t>(inj_name);
    EKAT_REQUIRE_MSG(
      tf.size() == nctrl,
      "Incorrect number of values in transfer function for injection site \"" + inj_name + "\""
    );
    m_inj_transfer_funcs.push_back(tf);

    // controller history
    // m_ctrl_hist.push_back(std::vector<double>());

  }

}

bool Controller::
is_ctrl_time(const int ymd, const int tod) {
  // ymd and tod should be for end of this time step, not previous
  // TODO: interval checking needs to be generalized (e.g., seasonal)

  TimeStamp curr_time(ymd, tod);

  // const int curr_year = curr_time.year();
  const int curr_month = curr_time.month();
  const int curr_day = curr_time.day();

  // controller calculations will only trigger at curr_day = 1 and tod == 0
  if ((curr_day == 1) && (tod == 0)) {
    // only within controller bounds
    if ((curr_time >= m_start_time) || (curr_time <= m_end_time)) {
      // curr_month is one indexed
      if (((curr_month - 1) % m_ctrl_interval_months) == 0) {
        return true;
      }
    }
  }

  return false;

}

void Controller::
operator()(const int ymd) {

  // check that controller hasn't already been calculated for this
  EKAT_REQUIRE_MSG(
    m_ctrl_hist.find(ymd) == m_ctrl_hist.end(),
    "Controller already computed for date: " + std::to_string(ymd)
  );

  // calculate controller
  std::vector<double> mass_inj = {1.0, 2.0};

  // record controller history
  m_ctrl_hist[ymd] = mass_inj;

}

std::vector<double> Controller::
get_ctrl_hist_at(const int ymd) { return m_ctrl_hist.at(ymd); }

}