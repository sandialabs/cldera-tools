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

} // namespace cldera

#endif // CLDERA_FIELD_STAT
