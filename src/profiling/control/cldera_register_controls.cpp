#include "cldera_control.hpp"

namespace cldera {

void register_controls ()
{
  auto& factory = ControlFactory::instance();
  //factory.register_product("scalar",&create_stat<ScalarControl>);
}

} // namespace cldera
