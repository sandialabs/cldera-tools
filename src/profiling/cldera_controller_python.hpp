#ifndef CLDERA_CONTROLLER_PYTHON_HPP
#define CLDERA_CONTROLLER_PYTHON_HPP

#include <vector>

namespace cldera {

class ProfilingContext;
class Controller;

void call_python_controller_hook(ProfilingContext& c,
                                 Controller& ctrl,
                                 const int ymd,
                                 const int tod,
                                 std::vector<double>& out_mass_inj);

void controller_finalize_python();

} // namespace cldera

#endif // CLDERA_CONTROLLER_PYTHON_HPP
