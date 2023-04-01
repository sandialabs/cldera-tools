#ifndef SRC_PROFILING_STATS_CLDERA_FIELD_ZONAL_MEAN_HPP_
#define SRC_PROFILING_STATS_CLDERA_FIELD_ZONAL_MEAN_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <memory>

namespace cldera {

class FieldZonalMean : public FieldStatAlongAxis
{
public:
  FieldZonalMean (const ekat::Comm& comm, const ekat::ParameterList& pl);

  std::string type () const override { return "zonal_mean"; }

  std::vector<std::string> get_aux_fields_names () const override {
    return {"lat", "area"};
  }

  void create_stat_field () override;

protected:

  void set_aux_fields_impl () override;

  void compute_impl () override;

  template <typename T, int N>
  void do_compute_impl ();

  const Bounds<Real> m_lat_bounds;
  const Bounds<int>  m_lev_bounds;

  Field m_lat, m_area;
  Real m_zonal_area = 0.0;

  std::vector<char> m_temp_memory;
};

} // namespace cldera

#endif /* SRC_PROFILING_STATS_CLDERA_FIELD_ZONAL_MEAN_HPP_ */
