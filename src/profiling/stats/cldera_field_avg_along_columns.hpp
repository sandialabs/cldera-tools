#ifndef CLDERA_FIELD_AVG_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_AVG_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_sum_along_columns.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class FieldAvgAlongColumns : public FieldSumAlongColumns
{
public:
  FieldAvgAlongColumns (const ekat::Comm& comm,
                        const ekat::ParameterList& pl)
   : FieldSumAlongColumns(comm,pl)
  { /* Nothing to do here */ }

  std::string type () const override { return "avg_along_columns"; }

protected:
  // NOTE: unlike global max/min/sum, we don't support IntType,
  //       so no need for extra template function
  void compute_impl () override {
    // Sum along columns
    FieldSumAlongColumns::compute_impl();

    // Divide by number of columns
    const auto& field_layout = m_field.layout();
    const auto& field_dims = field_layout.dims();
    const auto& field_names = field_layout.names();
    const auto it = std::find(field_names.begin(), field_names.end(), "ncol");
    long long size = field_dims[it-field_names.begin()];
    long long global_size;

    // Clock MPI ops
    track_mpi_all_reduce(m_comm,&size,&global_size,1,MPI_SUM, name());

    auto avg_field = m_stat_field.view_nonconst<Real>();
    int stat_size = avg_field.size();
    for (int i = 0; i < stat_size; ++i)
      avg_field(i) /= global_size;
  }
};

} // namespace cldera

#endif /* CLDERA_FIELD_AVG_ALONG_COLUMNS_HPP_ */
