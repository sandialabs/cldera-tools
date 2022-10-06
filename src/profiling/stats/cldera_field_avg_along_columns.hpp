#ifndef CLDERA_FIELD_AVG_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_AVG_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_sum_along_columns.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <limits>

namespace cldera {

class FieldAvgAlongColumns : public FieldSumAlongColumns
{
public:
  FieldAvgAlongColumns (const ekat::Comm& comm)
   : FieldSumAlongColumns(comm),
     m_comm (comm)
  { /* Nothing to do here */ }

  std::string name () const override { return "avg_along_columns"; }

protected:
  void compute_impl (const Field& f, Field& stat) const  override {
    // Sum along columns
    FieldSumAlongColumns::compute_impl(f,stat);

    // Divide by number of columns
    const auto& field_layout = f.layout();
    const auto& field_dims = field_layout.dims();
    const auto& field_names = field_layout.names();
    const auto it = std::find(field_names.begin(), field_names.end(), "ncol");
    long long size = field_dims[it-field_names.begin()];
    long long global_size;
    m_comm.all_reduce(&size, &global_size, 1, MPI_SUM);

    auto avg_field = stat.view_nonconst();
    auto stat_size = avg_field.size();
    for (int i = 0; i < stat_size; ++i)
      avg_field(i) /= global_size;
  }

  const ekat::Comm m_comm;
};

} // namespace cldera

#endif /* CLDERA_FIELD_AVG_ALONG_COLUMNS_HPP_ */
