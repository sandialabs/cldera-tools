
#include <iostream>

#include "cldera_controller.hpp"
#include "cldera_controller_python.hpp"
#include <iostream>

namespace cldera {

Controller::
Controller(const ekat::ParameterList& params)
{

  using vos_t = std::vector<std::string>;
  using vod_t = std::vector<double>;
  using voi_t = std::vector<int>;

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

  // get control interval definitions
  // NOTE: assumes a 365 day calendar
  auto intervals = ctrl_params.sublist("Intervals");
  for (auto it = intervals.params_names_cbegin(); it != intervals.params_names_cend(); ++it) {
    IntervalDef idef;
    idef.name = *it;
    auto interval_bounds = intervals.get<voi_t>(idef.name);
    idef.start_mmdd = interval_bounds[0];
    idef.end_mmdd   = interval_bounds[1];
    idef.span_days  = span_days_365_mmdd(idef.start_mmdd, idef.end_mmdd);
    m_intervals.push_back(idef);
  }

  // Hard requirement: intervals must not overlap (365-day calendar, inclusive bounds)
  auto mmdd_to_doy = [](const int mmdd) {
    const int mm = mmdd / 100;
    const int dd = mmdd % 100;
    return day_of_year_365(mm, dd);
  };
  auto overlaps = [&](const IntervalDef& a, const IntervalDef& b) {
    auto in_interval = [&](const IntervalDef& i, const int doy) {
      const int s = mmdd_to_doy(i.start_mmdd);
      const int e = mmdd_to_doy(i.end_mmdd);
      if (e >= s) {
        return doy >= s && doy <= e;
      }
      return doy >= s || doy <= e; // wraps year end
    };
    const int a_s = mmdd_to_doy(a.start_mmdd);
    const int a_e = mmdd_to_doy(a.end_mmdd);
    return in_interval(b, a_s) || in_interval(b, a_e) ||
           in_interval(a, mmdd_to_doy(b.start_mmdd)) ||
           in_interval(a, mmdd_to_doy(b.end_mmdd));
  };

  for (size_t i = 0; i < m_intervals.size(); ++i) {
    for (size_t j = i + 1; j < m_intervals.size(); ++j) {
      EKAT_REQUIRE_MSG(
        !overlaps(m_intervals[i], m_intervals[j]),
        "Intervals overlap: '" + m_intervals[i].name + "' and '" + m_intervals[j].name + "'."
      );
    }
  }

  // get controller definitions (grouped by interval)
  auto ctrl_intervals = ctrl_params.sublist("Controllers");
  for (auto it = ctrl_intervals.sublists_names_cbegin(); it != ctrl_intervals.sublists_names_cend(); ++it) {
    const std::string& interval_name = *it;
    EKAT_REQUIRE_MSG(
      std::any_of(m_intervals.begin(), m_intervals.end(),
                  [&](const IntervalDef& d){ return d.name == interval_name; }),
      "Controller interval \"" + interval_name + "\" not found in Intervals list"
    );

    auto interval_pl = ctrl_intervals.sublist(interval_name);

    for (auto jt = interval_pl.sublists_names_cbegin(); jt != interval_pl.sublists_names_cend(); ++jt) {

      const std::string& ctrl_name = *jt;
      auto ctrl_def = interval_pl.sublist(ctrl_name);
      ControllerDef cdef;
      cdef.name = ctrl_name;
      cdef.interval = interval_name;

      // type
      // TODO: this needs to be updated if using MPC
      auto ctrl_type = ctrl_def.get<std::string>("Type");
      EKAT_REQUIRE_MSG(
        ctrl_type == "constant" || ctrl_type == "P" || ctrl_type == "PI" || ctrl_type == "PD" || ctrl_type == "PID",
        "Invalid controller type for definition \"" + ctrl_name + "\" must be one of: P, PI, PD, PID"
      );
      cdef.type = ctrl_type;

      // simple constant controller for debugging
      if (ctrl_type  == "constant") {
        cdef.const_val = ctrl_def.get<double>("Value");

      // actual controller
      } else {

        // gains
        // TODO: need to be updated if using MPC
        auto ctrl_gains = ctrl_def.get<vod_t>("Gains");
        EKAT_REQUIRE_MSG(
          ctrl_gains.size() == ctrl_type.size(),
          "Incorrect number of gains for definition \"" + ctrl_name
        );
        cdef.gains = ctrl_gains;

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
        cdef.field = ctrl_field;
        cdef.stat  = ctrl_stat;

        // residual reference value
        auto ctrl_refval = ctrl_def.get<double>("Reference");
        cdef.refval = ctrl_refval;

        // TODO: some checking on projections (constant does not need a projection)
        auto ctrl_proj = ctrl_def.get<std::string>("Projection");
        cdef.proj = ctrl_proj;

        // TODO: this needs to be generalized greatly
        if (ctrl_def.isSublist("Feedforward")) {
          auto ff_def = ctrl_def.sublist("Feedforward");
          cdef.ff_type = ff_def.get<std::string>("Type");
          cdef.ff_val  = ff_def.get<double>("Value");
        }

        // constraints
        if (ctrl_def.isParameter("Constraints")) {
          auto cnstr_names = ctrl_def.get<vos_t>("Constraints");
          for (std::string& cnstr_name : cnstr_names) {
            auto cnstr_def = ctrl_def.get<std::string>("cnstr_" + cnstr_name);
            cdef.cnstr_defs.push_back(cnstr_def);
          }
        }

      } // end defining actual (non-constant) controller

      m_controllers.push_back(cdef);

    } // end per-controller

  } // end per-interval


  // injection definitions
  auto injections = ctrl_params.sublist("Injections");
  for (auto it = injections.sublists_names_cbegin(); it != injections.sublists_names_cend(); ++it) {
    const std::string& inj_name = *it;
    m_inj_names.push_back(inj_name);
    auto inj_def = injections.sublist(inj_name);

    // TODO: lat-lon-alt may not be used, with Graham's current mask file method
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

    // controller history
    // m_ctrl_hist.push_back(std::vector<double>());

  }

  // transfer function definitions
  // TODO: generalize to a transfer function for each control interval?
  if (ctrl_params.isSublist("Transfer Function")) {
    auto tf_defs = ctrl_params.sublist("Transfer Function");

    int tf_size;
    for (int inj_idx = 0; inj_idx < m_inj_names.size(); ++inj_idx) {
      std::string inj_name = m_inj_names[inj_idx];
      auto tf = tf_defs.get<vod_t>(inj_name);

      // make sure all transfer functions are the same size
      if (inj_idx == 0) {
        tf_size = tf.size();
      } else {
        EKAT_REQUIRE_MSG(
          tf.size() == tf_size,
          "Incorrect number of values in transfer function for injection site \"" + inj_name + "\""
        );
      }
      m_inj_transfer_funcs.push_back(tf);
    }

    // check that the transfer function size matches the number of controllers for all intervals
    {
      std::map<std::string,int> counts;
      for (const auto& c : m_controllers) {
        counts[c.interval] += 1;
      }
      EKAT_REQUIRE_MSG(!counts.empty(),
        "No controllers defined while processing transfer functions.");

      const int expected = counts.begin()->second;
      for (const auto& kv : counts) {
        EKAT_REQUIRE_MSG(kv.second == expected,
          "Intervals have different numbers of controllers. Interval '" +
          kv.first + "' has " + std::to_string(kv.second) +
          ", expected " + std::to_string(expected) + ".");
      }
      EKAT_REQUIRE_MSG(tf_size == expected,
        "Transfer function size (" + std::to_string(tf_size) +
        ") does not match number of controllers per interval (" +
        std::to_string(expected) + ").");
    }

  } else {
    EKAT_ERROR_MSG("Either no transfer function definition has been supplied, or the code has not been updated to permit per-interval transfer functions.");
  }
  

}

bool Controller::
is_ctrl_time(const int ymd, const int tod) {
  // ymd and tod should be for end of this time step, not previous
  // TODO: interval checking needs to be generalized (e.g., seasonal)

  TimeStamp curr_time(ymd, tod);

  const int curr_month = curr_time.month();
  const int curr_day = curr_time.day();
  const int curr_mmdd = curr_month * 100 + curr_day;

  // controller calculations will only trigger at the start of an interval at tod == 0
  if (tod == 0) {
    // only within controller bounds
    if ((curr_time >= m_start_time) && (curr_time <= m_end_time)) {
      for (const auto& idef : m_intervals) {
        if (curr_mmdd == idef.start_mmdd) {
          return true;
        }
      }
    }
  }

  return false;

}

std::vector<double> Controller::
operator()(const int ymd, const int tod) {

  // check that controller hasn't already been calculated for this
  EKAT_REQUIRE_MSG(
    m_ctrl_hist.find(ymd) == m_ctrl_hist.end(),
    "Controller already computed for date: " + std::to_string(ymd)
  );

  std::vector<double> mass_inj;
  if (m_context) {
    call_python_controller_hook(*m_context, *this, ymd, tod, mass_inj);
  }

  // record controller history
  m_ctrl_hist[ymd] = mass_inj;

  return mass_inj;
}

std::vector<double> Controller::
get_ctrl_hist_at(const int ymd) { return m_ctrl_hist.at(ymd); }

}
