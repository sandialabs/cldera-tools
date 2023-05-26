#ifndef CLDERA_FIELD_STAT_HPP
#define CLDERA_FIELD_STAT_HPP

#include "profiling/cldera_field.hpp"

#include "timing/cldera_timing_session.hpp"

#include <ekat/util/ekat_factory.hpp>
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class FieldStat
{
public:
  FieldStat (const ekat::Comm& comm,
             const ekat::ParameterList& pl)
   : m_params (pl)
   , m_comm (comm)
  {
    m_name = m_params.get("Name",pl.name());
  }

  virtual ~FieldStat () = default;

  // The name of this field stat
  std::string name () const { return m_name; }

  // Unlike the previous, this should be the same for all instances of the same type
  virtual std::string type () const = 0;

  // Given a field, return the layout that the computed stat will have
  virtual FieldLayout stat_layout (const FieldLayout& field_layout) const = 0;

  // If derived stats need auxiliary fields, they need to override these
  virtual std::vector<std::string> get_aux_fields_names () const { return {}; }
  void set_aux_fields (const std::map<std::string,Field>& fields);

  template<typename... Fs>
  void set_aux_fields (const Fs&... fields);

  void set_field (const Field& f);

  // Compute the stat field
  Field compute (const TimeStamp& timestamp);

  // NOTE: For most stats, the stat data type matches the field one, but it might not be.
  //       E.g., a stat that stores max location would have stat data type IntType,
  //       regardless of the field data type. So make method virtual, to allow flexibility.
  virtual DataType stat_data_type() const;

  const Field& get_stat_field () const { return m_stat_field; }

protected:
  // If derived classes need to perform some additional setup steps when setting the field,
  // they can override this method

  virtual void set_field_impl (const Field& /* f */) {}
  virtual void set_aux_fields_impl (const std::map<std::string,Field>&) {}
  virtual void compute_impl () = 0;

  ekat::ParameterList   m_params;
  ekat::Comm            m_comm;

  std::string           m_name;
  TimeStamp             m_timestamp;

  bool   m_aux_fields_set = false;
  Field  m_field;
  Field  m_stat_field;
};

template<typename... Fs>
void FieldStat::set_aux_fields (const Fs&... fields) {
  EKAT_REQUIRE_MSG ((ekat::SameType<Field,Fs...>::value),
      "Error! FieldStats::set_aux_fields needs a variadic list of of fields as input.\n");

  std::vector<Field> v { {fields...} };
  std::map<std::string,Field> aux_fs;
  for (const auto& f : v) {
    aux_fs[f.name()] = f;
  }
  set_aux_fields(aux_fs);
}

// ================= FACTORY for stat creation ================== //
using StatFactory =
  ekat::Factory<FieldStat,
                std::string,
                std::shared_ptr<FieldStat>,
                const ekat::Comm&,
                const ekat::ParameterList&>;

template<typename StatType>
inline std::shared_ptr<FieldStat>
create_stat (const ekat::Comm& comm, const ekat::ParameterList& params) {
  return std::make_shared<StatType>(comm,params);
}

// Special case of stat, returning a scalar
class FieldScalarStat : public FieldStat
{
public:
  FieldScalarStat (const ekat::Comm& comm,
                   const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  { /* Nothing to do here */ }

  FieldLayout stat_layout (const FieldLayout& /*field_layout*/) const override {
    return FieldLayout();
  }
};

// Special case of stat, returning a single part field
class FieldSinglePartStat : public FieldStat
{
public:
  FieldSinglePartStat (const ekat::Comm& comm,
                       const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  { /* Nothing to do here */ }

  std::vector<int> compute_stat_strides(const FieldLayout& field_layout) const;

  int compute_stat_index(const int ipart, const int part_index,
                         const int field_rank, const int field_part_dim,
                         const std::vector<int>& part_dims,
                         const std::vector<int>& stat_strides) const;
};

} // namespace cldera

#endif // CLDERA_FIELD_STAT
