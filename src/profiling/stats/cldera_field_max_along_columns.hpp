#ifndef CLDERA_FIELD_MAX_ALONG_COLUMNS_HPP_
#define CLDERA_FIELD_MAX_ALONG_COLUMNS_HPP_

#include "profiling/stats/cldera_field_stat_along_axis.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class FieldMaxAlongColumns : public FieldStatAlongAxis
{
public:
  FieldMaxAlongColumns (const ekat::Comm& comm,
                        const ekat::ParameterList& pl);

  std::string type () const override { return "max_along_columns"; }

protected:
  void compute_impl () override;

  template<typename T, int N>
  void do_compute_impl ();
};

} // namespace cldera

#endif /* CLDERA_FIELD_MAX_ALONG_COLUMNS_HPP_ */
