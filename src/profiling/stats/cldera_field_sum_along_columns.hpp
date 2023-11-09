#ifndef CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class FieldSumAlongColumns : public FieldStatAlongAxis
{
public:
  FieldSumAlongColumns (const ekat::Comm& comm,
                        const ekat::ParameterList& params);

  std::string type () const override { return "sum_along_columns"; }

  void create_stat_field () override;
protected:
  void compute_impl () override;

  template<typename T, int N>
  void do_compute_impl ();

  std::vector<char> m_scratch;
};

} // namespace cldera

#endif /* CLDERA_FIELD_SUM_ALONG_COLUMNS_HPP_ */
