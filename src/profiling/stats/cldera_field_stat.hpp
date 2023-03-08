#ifndef CLDERA_FIELD_STAT_HPP
#define CLDERA_FIELD_STAT_HPP

#include "profiling/cldera_field.hpp"

namespace cldera {

class FieldStat
{
public:
  virtual ~FieldStat () = default;

  // The name of this field stat
  virtual std::string name () const = 0;

  // Given a field, return the layout that the computed stat will have
  virtual FieldLayout stat_layout (const FieldLayout& field_layout) const = 0;

  // Compute the stat field
  void compute (const Field& f, Field& stat) const {

    // Sanity checks
    EKAT_REQUIRE_MSG (stat.layout()==stat_layout(f.layout()),
        "Error! Invalid layout for input stat field.\n");

    EKAT_REQUIRE_MSG (f.committed(), "Error! Input field is not committed.\n");
    EKAT_REQUIRE_MSG (stat.committed(), "Error! Input stat field is not committed.\n");

    compute_impl(f,stat);
  }

  // NOTE: For most stats, the stat data type matches the field one, but it might not be.
  //       E.g., a stat that stores max location would have stat data type IntType,
  //       regardless of the field data type. So make method virtual, to allow flexibility.
  virtual DataType stat_data_type(const Field& f) const {
    return f.data_type();
  }

  // Shortcut if you don't have a pre-built field
  Field compute (const Field& f) const {
    Field stat (f.name() + "_" + name(), stat_layout(f.layout()), DataAccess::Copy, stat_data_type(f));
    stat.commit();
    compute(f,stat);
    return stat;
  }

protected:
  virtual void compute_impl (const Field& f, Field& stat) const = 0;
};

// Special case of stat, returning a scalar
class FieldScalarStat : public FieldStat
{
public:
  FieldLayout stat_layout (const FieldLayout& /*field_layout*/) const override {
    return FieldLayout();
  }
};

// Special case of stat, returning a single part field
class FieldSinglePartStat : public FieldStat
{
public:
  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return FieldLayout(field_layout.dims(), field_layout.names());
  }

  inline std::vector<int> compute_stat_strides(const FieldLayout& field_layout) const
  {
    const int field_rank = field_layout.rank();
    const auto& field_dims = field_layout.dims();
    std::vector<int> stat_strides(field_rank, 1);
    for (int axis = field_rank-2; axis >= 0; --axis) { // layout right
      stat_strides[axis] = stat_strides[axis+1] * field_dims[axis+1];
    }
    return stat_strides;
  }

  inline int compute_stat_index(const int ipart, const int part_index, const int field_rank, const int field_part_dim,
      const std::vector<int>& part_dims, const std::vector<int>& stat_strides) const
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
};

} // namespace cldera

#endif // CLDERA_FIELD_STAT
