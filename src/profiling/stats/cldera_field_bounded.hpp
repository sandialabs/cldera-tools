#ifndef CLDERA_FIELD_BOUNDED_HPP_
#define CLDERA_FIELD_BOUNDED_HPP_

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldBounded : public FieldSinglePartStat
{
public:
  FieldBounded (const ekat::Comm& comm, const ekat::ParameterList& pl);

  std::string type () const override { return "bounded"; }

  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return field_layout;
  }

protected:
  void compute_impl () override;

  template <typename T, int N>
  void do_compute_impl ();

  const Bounds<Real> m_bounds;
  const Real m_mask_val;
};

} // namespace cldera

#endif /* CLDERA_FIELD_BOUNDED_HPP_ */
