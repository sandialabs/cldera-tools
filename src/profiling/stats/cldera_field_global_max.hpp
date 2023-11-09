#ifndef CLDERA_FIELD_GLOBAL_MAX_HPP
#define CLDERA_FIELD_GLOBAL_MAX_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>


namespace cldera {

class FieldGlobalMax : public FieldScalarStat
{
public:
  FieldGlobalMax (const ekat::Comm& comm,
                  const ekat::ParameterList& pl);

  std::string type () const override { return "global_max"; }

protected:
  void compute_impl () override;

  template<typename T, int N>
  void do_compute_impl ();
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_MAX_HPP
