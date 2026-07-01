#include "bindingsSymbolic.hpp"

#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "SEB.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"

py::object pyseb_expression_to_sympy(const Expression& expr)
{
    return py::module::import("pyseb.symbolic").attr("SymPyExpression")(expr.to_python());
}

py::object pyseb_zero_sympy_expression()
{
    return py::module::import("pyseb.symbolic").attr("SymPyExpression")("0");
}

void register_symbolic_world_bindings(py::class_<World>& world)
{
    world.def("FormFactor", [](World& self, const std::string& myself, int depth, int varForm) {
            return pyseb_expression_to_sympy(self.FormFactor(myself, depth, varForm));
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("FormFactor_Unnormalized", [](World& self, const std::string& myself, int depth, int varForm) {
            return pyseb_expression_to_sympy(self.FormFactor_Unnormalized(myself, depth, varForm));
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("FormFactorGeneric", [](World& self, const std::string& myself, int depth) {
            return pyseb_expression_to_sympy(self.FormFactorGeneric(myself, depth));
        }, py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH)
        .def("FormFactorAmplitude", [](World& self, refPoint& ref, int depth, int varForm) {
            return pyseb_expression_to_sympy(self.FormFactorAmplitude(ref, depth, varForm));
        }, py::arg("ref"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = static_cast<int>(QVAR))
        .def("PhaseFactor", [](World& self, refPoint& ref1, refPoint& ref2, int depth) {
            return pyseb_expression_to_sympy(self.PhaseFactor(ref1, ref2, depth));
        }, py::arg("ref1"), py::arg("ref2"), py::arg("depth") = WORLDMAXDEPTH)
        .def("RadiusOfGyration2", [](World& self, const std::string& myself) {
            return pyseb_expression_to_sympy(self.RadiusOfGyration2(myself));
        }, py::arg("myself"))
        .def("SMSD_ref2ref", [](World& self, refPoint& ref1, refPoint& ref2) {
            return pyseb_expression_to_sympy(self.SMSD_ref2ref(ref1, ref2));
        }, py::arg("ref1"), py::arg("ref2"))
        .def("SMSD_ref2scat", [](World& self, refPoint& ref) {
            return pyseb_expression_to_sympy(self.SMSD_ref2scat(ref));
        }, py::arg("ref"))
        .def("SMSD_scat2scat", [](World& self) {
            (void)self;
            return pyseb_zero_sympy_expression();
        });
}
