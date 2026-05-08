#ifndef CLDERA_FIELD_MASKED_WRITE_HPP
#define CLDERA_FIELD_MASKED_WRITE_HPP

#include "profiling/stats/cldera_field_stat.hpp"

namespace cldera {

class FieldMaskedWrite : public FieldStat {
public:

  FieldMaskedWrite (const ekat::Comm& comm,
                    const ekat::ParameterList& pl);

  std::string type () const override { return "masked_write"; }

  // Return a factor to convert Tg/yr to molecules/s (need to multiply by volume in cm^3 after this conversion)
  Real get_conversion_factor () const { return m_sec_per_year * m_g_per_tg / m_avogadros_number * m_molecular_weight; }

  // Set the write amounts according to a vector input - must match m_write_values size
  void set_write_values(const std::vector<Real> &write_values);

  std::vector<std::string> get_aux_fields_names () const override;

  FieldLayout stat_layout (const FieldLayout& fl) const override;

  // Writing scalar data to fields is real-valued
  DataType stat_data_type() const override { return DataType::RealType; }
protected:

  void set_aux_fields_impl () override;

  void compute_write_heights ();

  void compute_write_volumes ();

  void compute_impl () override;

  void print_vars ();

  template<typename T, int N>
  void do_compute_impl ();

  void load_mask_field (const Field& my_col_gids);

  // The mask field
  Field         m_mask_field;

  // The height field - zi
  Field         m_height_field;

  // The area field - area
  Field         m_area_field;

  // The number of injection sites
  int m_num_injection_sites;

  // Map every mask value to an index in [0, m_num_injection_sites)
  std::map<int,int>   m_mask_val_to_stat_entry;
  
  // m_mask_cols[imask][icol] gives the local column ids to write at the requested mask
  std::vector<std::vector<int>> m_mask_cols;

  // m_ilev[imask][icol] gives the correct ilev to write at the requested height
  std::vector<std::vector<int>> m_ilev;

  // m_mask_area[imask][icol] gives the area (not volume) of the mask column
  std::vector<std::vector<Real>> m_mask_area;

  // m_mask_thickness[imask][icol] gives the thickness of the mask column at the injection height
  std::vector<std::vector<Real>> m_mask_thickness;

  // m_mask_height[imask] gives the height of the write
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
  const Real m_avogadros_number = 6.02214076e23;
  const Real m_molecular_weight = 64.066;
  const Real m_sec_per_year = 31536000;
  const Real m_g_per_tg = 1.0e12;

  // Allows to have a stat that saves the mask field
  bool          m_output_mask_field;
};

} // namespace cldera

#endif // CLDERA_FIELD_MASKED_WRITE_HPP
