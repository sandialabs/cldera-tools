#ifndef CLDERA_FIELD_BOUNDED_MASKED_INTEGRAL_HPP
#define CLDERA_FIELD_BOUNDED_MASKED_INTEGRAL_HPP

#include "profiling/stats/cldera_field_masked_integral.hpp"

namespace cldera {

class FieldBoundedMaskedIntegral : public FieldMaskedIntegral {
public:

  FieldBoundedMaskedIntegral (const ekat::Comm& comm,
                              const ekat::ParameterList& pl);

  std::string type () const override { return "bounded_masked_integral"; }

protected:

  void compute_impl () override;

  template<typename T>
  void do_compute_impl ();

  // The bounds outside of which the field values must be discarded
  Bounds<Real>  m_bounds;
  bool          m_has_bounds;
};

} // namespace cldera

#endif // CLDERA_FIELD_BOUNDED_MASKED_INTEGRAL_HPP
