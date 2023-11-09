#ifndef CLDERA_FIELD_GLOBAL_MIN_HPP
#define CLDERA_FIELD_GLOBAL_MIN_HPP

#include "profiling/stats/cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class FieldGlobalMin : public FieldScalarStat
{
public:
  FieldGlobalMin (const ekat::Comm& comm,
                  const ekat::ParameterList& pl);

  std::string type () const override { return "global_min"; }

protected:
  void compute_impl () override;

  template<typename T, int N>
  void do_compute_impl ();
};

} // namespace cldera

#endif // CLDERA_FIELD_GLOBAL_MIN_HPP

