#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <sstream>

#include "SEB.hpp"
#include "World.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "SymbolicInterface.hpp"
#include "GiNaCSymbolic.hpp"

namespace py = pybind11;

// Helper function to convert GiNaC Expression to Python SymPyExpression
py::object convert_ginac_to_sympy(const Expression& expr) {
    return py::module::import("pyseb.symbolic").attr("SymPyExpression")(expr.to_string());
}

void register_ginac_bindings(py::class_<World>& world) {
    world.def("FormFactor", [](World& self, const std::string& myself, int depth, int varForm) {
            auto expr = self.FormFactor(myself, depth, varForm);
            return convert_ginac_to_sympy(expr);
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("FormFactor_Unnormalized", [](World& self, const std::string& myself, int depth, int varForm) {
            auto expr = self.FormFactor_Unnormalized(myself, depth, varForm);
            return convert_ginac_to_sympy(expr);
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("FormFactorGeneric", [](World& self, const std::string& myself, int depth) {
            auto expr = self.FormFactorGeneric(myself, depth);
            return convert_ginac_to_sympy(expr);
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH)
        .def("FormFactorAmplitude", [](World& self, refPoint& ref, int depth, int varForm) {
            auto expr = self.FormFactorAmplitude(ref, depth, varForm);
            return convert_ginac_to_sympy(expr);
        }, py::arg("ref"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("PhaseFactor", [](World& self, refPoint& ref1, refPoint& ref2, int depth) {
            auto expr = self.PhaseFactor(ref1, ref2, depth);
            return convert_ginac_to_sympy(expr);
        }, py::arg("ref1"), py::arg("ref2"), py::arg("depth") = WORLDMAXDEPTH)
        .def("RadiusOfGyration2", [](World& self, const std::string& myself) {
            auto expr = self.RadiusOfGyration2(myself);
            return convert_ginac_to_sympy(expr);
        }, py::arg("myself"))
        .def("SMSD_ref2ref", [](World& self, refPoint& ref1, refPoint& ref2) {
            auto expr = self.SMSD_ref2ref(ref1, ref2);
            return convert_ginac_to_sympy(expr);
        }, py::arg("ref1"), py::arg("ref2"))
        .def("SMSD_ref2scat", [](World& self, refPoint& ref) {
            auto expr = self.SMSD_ref2scat(ref);
            return convert_ginac_to_sympy(expr);
        }, py::arg("ref"))
        .def("SMSD_scat2scat", [](World& self) {
            return py::module::import("pyseb.symbolic").attr("SymPyExpression")("0");
        });
}
