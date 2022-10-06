#include "cldera_field_global_max.hpp"
#include "cldera_field_global_min.hpp"
#include "cldera_field_global_sum.hpp"
#include "cldera_field_global_avg.hpp"

#include <ekat/util/ekat_string_utils.hpp>

namespace cldera {

std::shared_ptr<FieldStat>
create_stat (const std::string& name, const ekat::Comm& comm) {
  std::shared_ptr<FieldStat> stat;
  ekat::CaseInsensitiveString name_ci = name;
  if (name_ci=="global_max") {
    stat = std::make_shared<FieldGlobalMax>(comm);
  } else if (name_ci=="global_min") {
    stat = std::make_shared<FieldGlobalMin>(comm);
  } else if (name_ci=="global_sum") {
    stat = std::make_shared<FieldGlobalSum>(comm);
  } else if (name_ci=="global_avg") {
    stat = std::make_shared<FieldGlobalAvg>(comm);
  } else {
    EKAT_ERROR_MSG ("Unrecognized/unsupported stat '" + name + "'.\n");
  }

  return stat;
}

} // namespace cldera
