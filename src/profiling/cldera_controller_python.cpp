#include "cldera_controller_python.hpp"

#include "cldera_controller.hpp"
#include "cldera_profiling_context.hpp"
#include "cldera_profiling_archive.hpp"
#include "stats/cldera_field_single_rank.hpp"
#include "cldera_field.hpp"
#include "cldera_field_layout.hpp"
#include "cldera_data_type.hpp"
#include "cldera_profiling_types.hpp"
#include "cldera_time_stamp.hpp"

#include <ekat/ekat_assert.hpp>
#include <ekat/ekat_parameter_list.hpp>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <map>
#include <vector>

namespace py = pybind11;

namespace cldera {

static py::scoped_interpreter* py_guard = nullptr;

PYBIND11_EMBEDDED_MODULE(sctrctrl_bindings, m) {
  py::enum_<cldera::DataType>(m, "DataType")
    .value("Invalid", cldera::DataType::Invalid)
    .value("RealType", cldera::DataType::RealType)
    .value("IntType", cldera::DataType::IntType)
    .export_values();

  py::class_<cldera::FieldLayout>(m, "FieldLayout")
    .def("dims", &cldera::FieldLayout::dims)
    .def("names", &cldera::FieldLayout::names)
    .def("rank", &cldera::FieldLayout::rank)
    .def("size", &cldera::FieldLayout::size)
    .def("to_string", &cldera::FieldLayout::to_string);

  py::class_<cldera::Field>(m, "Field")
    .def("name", &cldera::Field::name)
    .def("layout", &cldera::Field::layout, py::return_value_policy::reference_internal)
    .def("data_type", &cldera::Field::data_type)
    .def("nparts", &cldera::Field::nparts)
    .def("part_dim", &cldera::Field::part_dim)
    .def("part_offset", &cldera::Field::part_offset)
    .def("committed", &cldera::Field::committed);

  py::class_<cldera::Controller>(m, "Controller")
    .def("is_ctrl_time", &cldera::Controller::is_ctrl_time)
    .def("get_ctrl_hist_at", &cldera::Controller::get_ctrl_hist_at)
    .def("run", &cldera::Controller::operator());

  py::class_<cldera::FieldSingleRank, std::shared_ptr<cldera::FieldSingleRank>>(m, "FieldSingleRank")
    .def("name", &cldera::FieldSingleRank::name)
    .def("type", &cldera::FieldSingleRank::type)
    .def("get_stat_field", &cldera::FieldSingleRank::get_stat_field,
         py::return_value_policy::reference_internal)
    .def("stat_data_real_view", [](const cldera::FieldSingleRank& s) {
      const auto& f = s.get_stat_field();
      EKAT_REQUIRE_MSG(f.committed(),
          "[CLDERA] Stat field must be committed before accessing data.\n");
      EKAT_REQUIRE_MSG(f.data_type() == cldera::RealType,
          "[CLDERA] Stat field data_type must be RealType.\n");
      EKAT_REQUIRE_MSG(f.nparts() == 1,
          "[CLDERA] Stat field must have a single part.\n");
      auto view = f.view<const cldera::Real>();
      const ssize_t n = view.extent(0);
      auto base = py::cast(&s, py::return_value_policy::reference);
      auto arr = py::array_t<cldera::Real>(
        {n},
        {static_cast<ssize_t>(sizeof(cldera::Real))},
        view.data(),
        base
      );
      arr.attr("setflags")(false);
      return arr;
    })
    .def("stat_data_int_view", [](const cldera::FieldSingleRank& s) {
      const auto& f = s.get_stat_field();
      EKAT_REQUIRE_MSG(f.committed(),
          "[CLDERA] Stat field must be committed before accessing data.\n");
      EKAT_REQUIRE_MSG(f.data_type() == cldera::IntType,
          "[CLDERA] Stat field data_type must be IntType.\n");
      EKAT_REQUIRE_MSG(f.nparts() == 1,
          "[CLDERA] Stat field must have a single part.\n");
      auto view = f.view<const int>();
      const ssize_t n = view.extent(0);
      auto base = py::cast(&s, py::return_value_policy::reference);
      auto arr = py::array_t<int>(
        {n},
        {static_cast<ssize_t>(sizeof(int))},
        view.data(),
        base
      );
      arr.attr("setflags")(false);
      return arr;
    });
}

void call_python_controller_hook (ProfilingContext& c,
                                  Controller& ctrl,
                                  const int ymd,
                                  const int tod,
                                  std::vector<double>& out_mass_inj)
{

  // Non-root ranks don't need to call controller
  const auto& comm = c.get_comm();
  if (!comm.am_i_root()) {
    return;
  }

  static constexpr const char* module_name = "sctrctrl";
  static constexpr const char* function_name = "run_controller";

  // Spin up Python interpreter (lifetime managed via controller_finalize_python()).
  if (!py_guard) {
    py_guard = new py::scoped_interpreter();
  }

  try {
    py::gil_scoped_acquire gil;
    py::module_::import("sctrctrl_bindings");
    py::module_ mod = py::module_::import(module_name);
    py::object func = mod.attr(function_name);

    using stat_ptr_t = std::shared_ptr<FieldStat>;
    using requests_t = std::map<std::string,std::vector<stat_ptr_t>>;
    std::vector<std::shared_ptr<FieldSingleRank>> singlerank_stats;

    if (c.has_data("requests")) {
      auto& requests = c.get<requests_t>("requests");
      for (auto& it : requests) {
        for (auto& stat : it.second) {
          if (stat->type() != "singlerank") {
            continue;
          }

          // Explicitly downcast from FieldStat so that pybind can access FieldSingleRank methods
          auto sr = std::dynamic_pointer_cast<FieldSingleRank>(stat);
          if (sr) {
            singlerank_stats.push_back(sr);
          }
        }
      }
    }

    // Append cached static SingleRank stats (gids/area) if available
    if (c.has_data("controller_static_singlerank")) {
      auto& cached = c.get<std::vector<std::shared_ptr<FieldSingleRank>>>("controller_static_singlerank");
      for (auto& s : cached) {
        if (s) { singlerank_stats.push_back(s); }
      }
    }

    // Pull SingleRankStat pointers into a Python list
    py::list singlerank_list;
    for (auto& sr : singlerank_stats) {
      singlerank_list.append(py::cast(sr));
    }

    // Call controller interface, get back mass injections
    py::object result = func(py::cast(&ctrl, py::return_value_policy::reference),
                             singlerank_list, ymd, tod);
    if (!result.is_none()) {
      out_mass_inj = result.cast<std::vector<double>>();
    }
  } catch (const py::error_already_set& e) {
    EKAT_ERROR_MSG(std::string("[CLDERA] Python error in controller hook:\n") + e.what());
  }
}

void controller_finalize_python()
{
  if (py_guard) {
    delete py_guard;
    py_guard = nullptr;
  }
}

} // namespace cldera
