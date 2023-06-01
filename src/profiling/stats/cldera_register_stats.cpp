#include "cldera_field_global_max.hpp"
#include "cldera_field_global_min.hpp"
#include "cldera_field_global_sum.hpp"
#include "cldera_field_global_avg.hpp"
#include "cldera_field_identity.hpp"
#include "cldera_field_max_along_columns.hpp"
#include "cldera_field_min_along_columns.hpp"
#include "cldera_field_sum_along_columns.hpp"
#include "cldera_field_avg_along_columns.hpp"
#include "cldera_field_bounded.hpp"
#include "cldera_field_bounding_box.hpp"
#include "cldera_field_pnetcdf_reference.hpp"
#include "cldera_field_zonal_mean.hpp"
#include "cldera_field_vertical_contraction.hpp"
#include "cldera_field_stat_pipe.hpp"
#include "cldera_field_masked_integral.hpp"

namespace cldera {

void register_stats ()
{
  auto& factory = StatFactory::instance();
  factory.register_product("global_max",&create_stat<FieldGlobalMax>);
  factory.register_product("global_min",&create_stat<FieldGlobalMin>);
  factory.register_product("global_sum",&create_stat<FieldGlobalSum>);
  factory.register_product("global_avg",&create_stat<FieldGlobalAvg>);

  factory.register_product("max_along_columns",&create_stat<FieldMaxAlongColumns>);
  factory.register_product("min_along_columns",&create_stat<FieldMinAlongColumns>);
  factory.register_product("sum_along_columns",&create_stat<FieldSumAlongColumns>);
  factory.register_product("avg_along_columns",&create_stat<FieldAvgAlongColumns>);

  factory.register_product("identity",&create_stat<FieldIdentity>);

  factory.register_product("bounded",&create_stat<FieldBounded>);
  factory.register_product("bounding_box",&create_stat<FieldBoundingBox>);
  factory.register_product("zonal_mean",&create_stat<FieldZonalMean>);

  factory.register_product("pnetcdf_reference",&create_stat<FieldPnetcdfReference>);
  factory.register_product("vertical_contraction",&create_stat<FieldVerticalContraction>);
  factory.register_product("pipe",&create_stat<FieldStatPipe>);
  factory.register_product("masked_integral",&create_stat<FieldMaskedIntegral>);
}

} // namespace cldera
