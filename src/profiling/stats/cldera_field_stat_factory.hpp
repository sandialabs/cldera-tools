#ifndef CLDERA_STATS_FACTORY_HPP
#define CLDERA_STATS_FACTORY_HPP

#include "cldera_field_stat.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <memory>
#include <string>

namespace cldera {

std::shared_ptr<FieldStat>
create_stat (const std::string& name, const ekat::Comm& comm);

} // namespace cldera

#endif // CLDERA_STATS_FACTORY_HPP
