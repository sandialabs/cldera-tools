#ifndef CLDERA_FIELD_GLOBAL_SUM_HPP
#define CLDERA_FIELD_GLOBAL_SUM_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>
#include <ekat/ekat_parameter_list.hpp>

namespace cldera {

class FieldGlobalSum : public FieldScalarStat
{
public:
  FieldGlobalSum (const ekat::Comm& comm,
                  const ekat::ParameterList& pl);

  std::string type () const override { return "global_sum"; }

protected:
  void compute_impl () override;
  
  template<typename T, int N>
  void do_compute_impl ();
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_SUM_HPP
