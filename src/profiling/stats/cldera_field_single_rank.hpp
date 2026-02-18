#ifndef CLDERA_FIELD_SINGLE_RANK_HPP
#define CLDERA_FIELD_SINGLE_RANK_HPP

#include "cldera_field_stat.hpp"

namespace cldera
{

class FieldSingleRank : public FieldStat
{
public:
  FieldSingleRank (const ekat::Comm& comm,
                 const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  { /* Nothing to do here */ }

  std::string type () const override { return "singlerank"; }

  // Given a field, return the layout that the computed stat will have
  void create_stat_field ();

  // Given a field, return the layout that the computed stat will have
  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return field_layout;
  }

  void compute_impl () override;
private:

  template<typename T,int N>
  void do_compute_impl ();

  std::vector<int> m_part_sizes;
  int m_local_size;
  int m_global_size;
};

} // namespace cldera

#endif // CLDERA_FIELD_SINGLE_RANK_HPP
