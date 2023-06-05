#include "profiling/stats/cldera_field_stat.hpp"

#include "timing/cldera_timing_session.hpp"

namespace cldera {

void FieldStat::
set_aux_fields (const std::map<std::string,Field>& fields) {
  EKAT_REQUIRE_MSG (m_field.committed(),
      "Error! Aux fields should be set *after* the input field.\n"
      " - stat name: " + name()  + "\n");

  for (const auto& it : get_aux_fields_names()) {
    EKAT_REQUIRE_MSG (fields.count(it)==1,
        "Error! Input map is missing a required auxiliary field:\n"
        " - stat name: " + name() + "\n"
        " - aux field name: " + it + "\n");
    EKAT_REQUIRE_MSG (fields.at(it).committed(),
        "Error! Input aux field is not yet committed.\n"
        " - stat name: " + name() + "\n"
        " - aux field name: " + it + "\n");
  }
  set_aux_fields_impl (fields);

  m_aux_fields_set = true;
}

void FieldStat::
set_field (const Field& f) {
  EKAT_REQUIRE_MSG (f.committed(),
      "Error! Field must be committed before being set in the stat.\n"
      " - stat name : " + name() + "\n"
      " - field name: " + f.name() + "\n");
  m_field = f;

  set_field_impl (f);
}

void FieldStat::
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
Field FieldStat::
compute (const TimeStamp& timestamp) {
  auto& ts = timing::TimingSession::instance();
  ts.start_timer ("profiling::compute_" + name());
  EKAT_REQUIRE_MSG (m_stat_field.committed(),
      "Error! Field must be set in the stat before calling compute.\n"
      " - stat name : " + name() + "\n");

  // Store timestamp, in case it's needed by derived class
  m_timestamp = timestamp;

  // Call derived class impl
  compute_impl();

  ts.stop_timer ("profiling::compute_" + name());
  return m_stat_field;
}

DataType FieldStat::
stat_data_type() const {
  EKAT_REQUIRE_MSG (m_field.committed(),
      "Error! Field must be set in the stat before calling stat_data_type.\n"
      " - stat name : " + name() + "\n");
  return m_field.data_type();
}

std::vector<int>
FieldSinglePartStat::
compute_stat_strides(const FieldLayout& field_layout) const
{
  const int field_rank = field_layout.rank();
  const auto& field_dims = field_layout.dims();
  std::vector<int> stat_strides(field_rank, 1);
  for (int axis = field_rank-2; axis >= 0; --axis) { // layout right
    stat_strides[axis] = stat_strides[axis+1] * field_dims[axis+1];
  }
  return stat_strides;
}

int FieldSinglePartStat::
compute_stat_index(const int ipart, const int part_index,
                   const int field_rank, const int field_part_dim,
                   const std::vector<int>& part_dims,
                   const std::vector<int>& stat_strides) const
{
  int field_ijk, work_index = part_index, stat_index = 0;
  for (int axis = field_rank-1; axis >= 0; --axis) { // layout right
    if (axis == field_part_dim) {
      int part_size = part_dims[field_part_dim];
      field_ijk = ipart * part_size + work_index % part_size;
    }
    else
      field_ijk = work_index % part_dims[axis];
    stat_index += field_ijk * stat_strides[axis];
    work_index /= part_dims[axis];
  }
  return stat_index;
}

} // namespace cldera
