#include "profiling/control/cldera_control.hpp"
#include "profiling/stats/cldera_field_stat.hpp"

#include "timing/cldera_timing_session.hpp"

namespace cldera {

// ------------------ Control implementation ------------------ //

void Control::
set_aux_fields (const std::map<std::string,Field>& fields) {
  EKAT_REQUIRE_MSG (m_field.committed(),
      "Error! Aux fields should be set *after* the input field.\n"
      " - stat name: " + name()  + "\n");

  for (const auto& it : get_aux_fields_names()) {
    if (fields.count(it)>0) {
      EKAT_REQUIRE_MSG (fields.at(it).committed(),
          "Error! Input aux field is not yet committed.\n"
          " - stat name: " + name() + "\n"
          " - aux field name: " + it + "\n");

      m_aux_fields[it] = fields.at(it).read_only();
    }
  }
  set_aux_fields_impl ();

  m_aux_fields_set = true;
}

void Control::
set_field (const Field& f) {
  EKAT_REQUIRE_MSG (f.committed(),
      "Error! Field must be committed before being set in the stat.\n"
      " - stat name : " + name() + "\n"
      " - field name: " + f.name() + "\n");
  m_field = f;

  set_field_impl (f);
}

void Control::
create_stat_field () {
  EKAT_REQUIRE_MSG (m_field.committed(),
      "Error! Cannot create stat field until input field is set.\n"
      " - stat name: " + name() + "\n");
  EKAT_REQUIRE_MSG (m_aux_fields_set or get_aux_fields_names().size()==0,
      "Error! Cannot create stat field until all aux fields are set.\n"
      " - stat name: " + name() + "\n");
  m_stat_field = Field (name(), stat_layout(m_field.layout()), DataAccess::Copy, stat_data_type());
  m_stat_field.commit();
}

// Compute the stat field
Field Control::
compute (const TimeStamp& timestamp) {
  auto& ts = timing::TimingSession::instance();
  ts.start_timer ("control::compute_" + name());
  EKAT_REQUIRE_MSG (m_stat_field.committed(),
      "Error! Field must be set in the stat before calling compute.\n"
      " - stat name : " + name() + "\n");

  // Store timestamp, in case it's needed by derived class
  m_timestamp = timestamp;

  // Call derived class impl
  compute_impl();

  ts.stop_timer ("control::compute_" + name());
  return m_stat_field;
}

DataType Control::
stat_data_type() const {
  EKAT_REQUIRE_MSG (m_field.committed(),
      "Error! Field must be set in the stat before calling stat_data_type.\n"
      " - stat name : " + name() + "\n");
  return m_field.data_type();
}

void Control::
check_aux_fields (const std::vector<std::string>& names) const
{
  auto map_keys = [&]() {
    std::string s;
    if (m_aux_fields.size()==0) {
      return s;
    }
    auto it = m_aux_fields.begin();
    s += it->first;
    for (it++; it!=m_aux_fields.end(); ++it) {
      s += ", " + it->first;
    }
    return s;
  };
  for (const auto& n : names) {
    EKAT_REQUIRE_MSG (m_aux_fields.count(n)>0,
        "Error! A required field was not set in the aux field map.\n"
        "  - missing field: " + n + "\n"
        "  - aux fields   : " + map_keys() + "\n");
  }
}

} // namespace cldera
