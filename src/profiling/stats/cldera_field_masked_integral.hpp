#ifndef CLDERA_FIELD_MASKED_INTEGRAL_HPP
#define CLDERA_FIELD_MASKED_INTEGRAL_HPP

#include "profiling/stats/cldera_field_stat.hpp"

namespace cldera {

class FieldMaskedIntegral : public FieldStat {
public:

  FieldMaskedIntegral (const ekat::Comm& comm,
                       const ekat::ParameterList& pl);

  std::string type () const { return "masked_integral"; }

  std::vector<std::string> get_aux_fields_names () const override;

  FieldLayout stat_layout (const FieldLayout& fl) const override;

  // Since we may have weights, let's just always use Real for the result.
  DataType stat_data_type() const { return DataType::RealType; }
protected:

  void set_aux_fields_impl (const std::map<std::string,Field>& fields) override;

  void compute_impl () override;

  template<typename T>
  void do_compute_impl ();

  void load_mask_field (const Field& my_col_gids);

  // The mask field
  Field         m_mask_field;

  // Map every mask value to an index in [0,N), with N=number_of_mask_values
  std::map<int,int>   m_mask_val_to_stat_entry;
  
  // Optionally, we weigh the integrand by a weight field
  bool          m_use_weight;
  bool          m_average;
  Field         m_weight_field;
  Field         m_weight_integral;
};

} // namespace cldera

#endif // CLDERA_FIELD_MASKED_INTEGRAL_HPP
