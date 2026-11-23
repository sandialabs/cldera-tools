
#include <algorithm>
#include <cmath>
#include <iostream>

#include "cldera_controller.hpp"
#include "cldera_controller_python.hpp"

namespace cldera {

namespace {
// BMW: the yaml parser stores "1" as int and "1.0" as double. Accept both for real params,
//      so that e.g. 'Value: 0' does not throw a cryptic any_cast error.
double get_real (const ekat::ParameterList& pl, const std::string& name) {
  EKAT_REQUIRE_MSG (pl.isParameter(name),
      "[Controller] Error! Parameter '" + name + "' not found in '" + pl.name() + "'.\n");
  if (pl.isType<int>(name)) {
    return static_cast<double>(pl.get<int>(name));
  }
  return pl.get<double>(name);
}
} // anonymous namespace

Controller::
Controller(const ekat::ParameterList& params)
{
  init(params);
}

Controller::
Controller(const ekat::ParameterList& params, const TimeStamp& case_t0)
 : m_case_t0(case_t0)
 , m_has_case_t0(true)
{
  init(params);
}

std::vector<std::string> Controller::
get_feedback_stats() const
{
  vos_t names;
  for (const auto& c : m_controllers) {
    // Only P/PI/... controllers use feedback. Constant controllers may name a stat just
    // to have it logged; that does not need annual means.
    if (c.type=="constant" or c.stat.empty()) continue;
    if (std::find(names.begin(),names.end(),c.stat)==names.end()) {
      names.push_back(c.stat);
    }
  }
  return names;
}

void Controller::
init(const ekat::ParameterList& params)
{

  using vos_t = std::vector<std::string>;
  using vod_t = std::vector<double>;
  using voi_t = std::vector<int>;

	// load controller parameters
  const auto prof_fields = params.get<vos_t>("Fields To Track");
	const auto ctrl_params_filename = params.get<std::string>("Controller");
  m_filename = ctrl_params_filename;
  auto ctrl_params = ekat::parse_yaml_file(ctrl_params_filename);

	// get time stamps
  // BMW: the start can be given either as an absolute date ('Start Day', 'Start Time'),
  //      or relative to the case start ('Start Years After Case Start: N' -> Jan 1 00:00
  //      of case start year + N). The relative form lets the same yaml be used for any
  //      case start date (e.g. "start the controller 1 full year after the case starts").
  if (ctrl_params.isParameter("Start Years After Case Start")) {
    EKAT_REQUIRE_MSG (not ctrl_params.isParameter("Start Day"),
        "[Controller] Error! Specify only one of 'Start Day' and 'Start Years After Case Start'.\n");
    EKAT_REQUIRE_MSG (m_has_case_t0,
        "[Controller] Error! 'Start Years After Case Start' requires the case start time.\n");
    m_start_years_after = ctrl_params.get<int>("Start Years After Case Start");
    EKAT_REQUIRE_MSG (m_start_years_after>=0,
        "[Controller] Error! 'Start Years After Case Start' must be >= 0.\n");
    // If the case does not start on Jan 1 at 00:00, the first full year starts the next Jan 1.
    int y0 = m_case_t0.year();
    if (not (m_case_t0.month()==1 and m_case_t0.day()==1 and m_case_t0.tod()==0)) {
      ++y0;
    }
    m_start_time = TimeStamp((y0+m_start_years_after)*10000 + 101, 0);
    m_start_mode = "Start Years After Case Start";
  } else {
    const int start_ymd = ctrl_params.get<int>("Start Day");
    const int start_tod = ctrl_params.get<int>("Start Time");
    m_start_time = TimeStamp(start_ymd, start_tod);
    m_start_mode = "Start Day";
  }
	const int end_ymd = ctrl_params.get<int>("End Day");
	const int end_tod = ctrl_params.get<int>("End Time");
	m_end_time = TimeStamp(end_ymd, end_tod);
  EKAT_REQUIRE_MSG (m_start_time<=m_end_time,
      "[Controller] Error! Controller start (" + m_start_time.to_string() +
      ") is after its end (" + m_end_time.to_string() + ").\n");

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
        cdef.const_val = get_real(ctrl_def,"Value");

        // BMW: a constant controller may optionally name a feedback Field/Stat. It does not
        //      use it to compute its output, but it is passed to (and logged by) the Python
        //      controller, which lets us test the full feedback data path before switching
        //      to a P/PI/PID controller.
        if (ctrl_def.isParameter("Stat")) {
          auto ctrl_field = ctrl_def.get<std::string>("Field");
          EKAT_REQUIRE_MSG(
            std::find(prof_fields.begin(), prof_fields.end(), ctrl_field) != prof_fields.end(),
            "Controller field \"" + ctrl_field + "\" not found in tracked fields"
          );
          auto field_stats = params.sublist(ctrl_field).get<vos_t>("Compute Stats");
          auto ctrl_stat = ctrl_def.get<std::string>("Stat");
          EKAT_REQUIRE_MSG(
            std::find(field_stats.begin(), field_stats.end(), ctrl_stat) != field_stats.end(),
            "Invalid controller statistic \"" + ctrl_stat + "\" for tracked field \"" + ctrl_field + "\"\n"
          );
          cdef.field = ctrl_field;
          cdef.stat  = ctrl_stat;
        }

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
        // BMW: optional. If missing, the Python controller computes it from the annual
        //      mean of the year before the controller starts, plus 'Reference Offset'.
        if (ctrl_def.isParameter("Reference")) {
          cdef.refval = get_real(ctrl_def,"Reference");
          cdef.has_refval = true;
        }
        if (ctrl_def.isParameter("Reference Offset")) {
          cdef.ref_offset = get_real(ctrl_def,"Reference Offset");
        }
        if (ctrl_def.isParameter("Initial Annual Mean")) {
          cdef.init_annual_mean = get_real(ctrl_def,"Initial Annual Mean");
        }
        // BMW: feedforward-only years before 'Feedback Start Year'
        if (ctrl_def.isParameter("Feedback Start Year")) {
          cdef.fb_start_year = ctrl_def.get<int>("Feedback Start Year");
        }
        if (ctrl_def.isParameter("Minimum Output")) {
          cdef.min_out = get_real(ctrl_def,"Minimum Output");
        }
        if (ctrl_def.isParameter("Maximum Output")) {
          cdef.max_out = get_real(ctrl_def,"Maximum Output");
        }
        if (ctrl_def.isParameter("Minimum Annual Coverage")) {
          cdef.min_coverage = get_real(ctrl_def,"Minimum Annual Coverage");
          EKAT_REQUIRE_MSG (cdef.min_coverage>0 and cdef.min_coverage<=1,
              "[Controller] Error! 'Minimum Annual Coverage' must be in (0,1].\n");
        }

        // TODO: some checking on projections (constant does not need a projection)
        auto ctrl_proj = ctrl_def.get<std::string>("Projection");
        cdef.proj = ctrl_proj;

        // TODO: this needs to be generalized greatly
        if (ctrl_def.isSublist("Feedforward")) {
          auto& ff_def = ctrl_def.sublist("Feedforward");
          cdef.ff_type = ff_def.get<std::string>("Type");
          cdef.ff_val  = get_real(ff_def,"Value");
          EKAT_REQUIRE_MSG (cdef.ff_type=="linear" or cdef.ff_type=="constant",
              "[Controller] Error! Feedforward 'Type' must be 'linear' or 'constant'.\n"
              " - controller: " + ctrl_name + "\n - type: " + cdef.ff_type + "\n");
          if (ff_def.isParameter("Start Year")) {
            cdef.ff_start_year = ff_def.get<int>("Start Year");
          }
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

bool Controller::
is_interval_start(const int ymd, const int tod) const {
  // BMW: true at tod==0 on the start day of ANY interval, regardless of the controller
  //      active window. The cldera interface triggers the (Python) controller at these
  //      times, and the Python controller returns zeros outside the active window. This
  //      way injection stops at 'End Day', and stays zero before the start (it would
  //      otherwise keep the last rates until the next trigger).
  if (tod != 0) return false;
  const int curr_mmdd = ymd % 10000;
  for (const auto& idef : m_intervals) {
    if (curr_mmdd == idef.start_mmdd) {
      return true;
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
