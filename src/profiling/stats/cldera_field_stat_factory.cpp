#include "cldera_field_global_max.hpp"
#include "cldera_field_global_min.hpp"
#include "cldera_field_global_sum.hpp"
#include "cldera_field_global_avg.hpp"
#include "cldera_field_identity.hpp"
#include "cldera_field_max_along_columns.hpp"
#include "cldera_field_min_along_columns.hpp"
#include "cldera_field_sum_along_columns.hpp"
#include "cldera_field_avg_along_columns.hpp"

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
  } else if (name_ci=="max_along_columns") {
    stat = std::make_shared<FieldMaxAlongColumns>(comm);
  } else if (name_ci=="min_along_columns") {
    stat = std::make_shared<FieldMinAlongColumns>(comm);
  } else if (name_ci=="sum_along_columns") {
    stat = std::make_shared<FieldSumAlongColumns>(comm);
  } else if (name_ci=="avg_along_columns") {
    stat = std::make_shared<FieldAvgAlongColumns>(comm);
  } else if (name_ci=="identity") {
    stat = std::make_shared<FieldIdentity>();
  } else {
    EKAT_ERROR_MSG ("Unrecognized/unsupported stat '" + name + "'.\n");
  }

  return stat;
}

} // namespace cldera
