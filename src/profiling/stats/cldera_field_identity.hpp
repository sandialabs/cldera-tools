#ifndef CLDERA_FIELD_IDENTITY_HPP
#define CLDERA_FIELD_IDENTITY_HPP

#include "cldera_field_stat.hpp"

namespace cldera
{

class FieldIdentity : public FieldStat
{
public:
  FieldIdentity (const ekat::Comm& comm,
                 const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  { /* Nothing to do here */ }

  std::string type () const override { return "identity"; }

  // Given a field, return the layout that the computed stat will have
  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return field_layout;
  }

  void compute_impl () override;
private:

  template<typename T,int N>
  void do_compute_impl ();
};

} // namespace cldera

#endif // CLDERA_FIELD_IDENTITY_HPP
