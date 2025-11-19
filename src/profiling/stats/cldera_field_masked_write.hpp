#ifndef CLDERA_FIELD_MASKED_WRITE_HPP
#define CLDERA_FIELD_MASKED_WRITE_HPP

#include "profiling/stats/cldera_field_stat.hpp"

namespace cldera {

class FieldMaskedWrite : public FieldStat {
public:

  FieldMaskedWrite (const ekat::Comm& comm,
                    const ekat::ParameterList& pl);

  std::string type () const override { return "masked_write"; }

  // Return a factor to convert Tg to molecules
  Real get_conversion_factor () const { return 1e9 * m_avogadro_number / m_molecular_weight; }

  std::vector<std::string> get_aux_fields_names () const override;

  FieldLayout stat_layout (const FieldLayout& fl) const override;

  // Writing scalar data to fields is real-valued
  DataType stat_data_type() const override { return DataType::RealType; }
protected:

  void set_aux_fields_impl () override;

  void compute_impl () override;

  template<typename T, int N>
  void do_compute_impl ();

  void load_mask_field (const Field& my_col_gids);

  // The mask field
  Field         m_mask_field;

  // The height field - zi
  Field         m_height_field;

  // The area field - area
  Field         m_area_field;

  // Map every mask value to an index in [0,N), with N=number_of_mask_values
  std::map<int,int>   m_mask_val_to_stat_entry;
  
  // mask_cols[imask][icol] gives the local column ids to write at the requested mask
  std::vector<std::vector<int>> m_mask_cols;

  // ilev[imask][icol] gives the correct ilev to write at the requested height
  std::vector<std::vector<int>> m_ilev;

  // mask_area[imask][icol] gives the area (not volume) of the mask column
  std::vector<Real> m_mask_area;

  //
  std::vector<Real> m_mask_height;

  /// Values obtained from the input deck
  // An array of scalars to write (in Tg/yr) at each of the nonzero-indexed write locations
  std::vector<Real> m_write_values;

  // An array of scalars defining the height to write to (in km)
  std::vector<Real> m_write_heights;

  // The default value to write at all zero-indexed write locations
  Real m_default_write;

  // The name of the mask field in the mask file
  std::string m_mask_field_name;

  // Constants for conversions
  const Real m_avogadro_number = 6.02214076e23; 
  const Real m_molecular_weight = 64 * 1e-3; // technically 64.066 - not sure this matters, but it's more sig figs

  // Allows to have a stat that saves the mask field
  bool          m_output_mask_field;
};

} // namespace cldera

#endif // CLDERA_FIELD_MASKED_WRITE_HPP
