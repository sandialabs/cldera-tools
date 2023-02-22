#ifndef CLDERA_STATS_FACTORY_HPP
#define CLDERA_STATS_FACTORY_HPP

#include "cldera_field_stat.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <memory>

namespace cldera {

std::shared_ptr<FieldStat>
create_stat (const ekat::ParameterList& pl, const ekat::Comm& comm);

} // namespace cldera

#endif // CLDERA_STATS_FACTORY_HPP
