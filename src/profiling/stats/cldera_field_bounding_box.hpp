#ifndef CLDERA_FIELD_BOUNDING_BOX_HPP
#define CLDERA_FIELD_BOUNDING_BOX_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class FieldBoundingBox : public FieldSinglePartStat
{
public:
  FieldBoundingBox (const ekat::Comm& comm,
                    const ekat::ParameterList& pl);

  std::string type () const override { return "bounding_box"; }

  FieldLayout stat_layout (const FieldLayout& field_layout) const override {
    return field_layout;
  }

  std::vector<std::string> get_aux_fields_names () const override {
    return {"lat", "lon"};
  }
protected:
  void set_aux_fields_impl () override;

  void compute_impl () override;

  template <typename T, int N>
  void do_compute_impl ();

  const Bounds<Real> m_lat_bounds, m_lon_bounds, m_lev_bounds;
  const Real m_mask_val;
  const ekat::Comm m_comm;
  Field m_lat, m_lon;
};

} // namespace cldera

#endif // CLDERA_FIELD_BOUNDING_BOX_HPP
