#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "SEB.hpp"
#include "World.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "SymbolInterface.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pyseb, m) {
    m.doc() = "Python bindings for the Scattering Equation Builder (SEB) library";
    
    // Expose GraphID type
    py::class_<GraphID>(m, "GraphID")
        .def(py::init<>())
        .def("__repr__", [](const GraphID &id) {
            return "GraphID(" + std::to_string(id) + ")";
        });
    
    // Expose World class
    py::class_<World>(m, "World")
        .def(py::init<std::string>(), py::arg("id") = "World")
        .def("Add", py::overload_cast<SubUnit*, subName, std::string>(&World::Add),
             py::arg("sub"), py::arg("name"), py::arg("tag") = "")
        .def("Add", py::overload_cast<std::string, subName, std::string>(&World::Add),
             py::arg("subunit_type"), py::arg("name"), py::arg("tag") = "")
        .def("Add", py::overload_cast<GraphID, structName>(&World::Add),
             py::arg("gid"), py::arg("name"))
        .def("Link", py::overload_cast<SubUnit*, refPoint, refPoint, std::string>(&World::Link),
             py::arg("sub"), py::arg("newr"), py::arg("oldr"), py::arg("tag") = "")
        .def("Link", py::overload_cast<std::string, refPoint, refPoint, std::string>(&World::Link),
             py::arg("subunit_type"), py::arg("newr"), py::arg("oldr"), py::arg("tag") = "")
        .def("Link", py::overload_cast<GraphID, refPoint, refPoint>(&World::Link),
             py::arg("gid"), py::arg("r1"), py::arg("r2"))
        .def("FormFactor", &World::FormFactor,
             py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = FORM)
        .def("FormFactor_Unnormalized", &World::FormFactor_Unnormalized,
             py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = FORM)
        .def("FormFactorGeneric", &World::FormFactorGeneric,
             py::arg("myself"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = FORM)
        .def("FormFactorAmplitude", &World::FormFactorAmplitude,
             py::arg("ref"), py::arg("depth") = WORLDMAXDEPTH, py::arg("varForm") = FORM)
        .def("PhaseFactor", &World::PhaseFactor,
             py::arg("ref1"), py::arg("ref2"), py::arg("depth") = WORLDMAXDEPTH)
        .def("getParams", &World::getParams)
        .def("getParamsq", &World::getParamsq)
        .def("getq", &World::getq, py::arg("value") = 0)
        .def("q", &World::q, py::arg("value") = 0)
        .def("setParameter", &World::setParameter,
             py::arg("pl"), py::arg("str"), py::arg("x"))
        .def("Count", &World::Count,
             py::arg("ref"), py::arg("depth") = WORLDMAXDEPTH)
        .def("CountPairs", &World::CountPairs,
             py::arg("name"), py::arg("depth") = WORLDMAXDEPTH)
        .def("Evaluate", py::overload_cast<ex, ParameterList&>(&World::Evaluate),
             py::arg("e"), py::arg("pl"))
        .def("Evaluate", py::overload_cast<ex, ParameterList&, double>(&World::Evaluate),
             py::arg("e"), py::arg("pl"), py::arg("q"))
        .def("Evaluate", py::overload_cast<ex, ParameterList&, DoubleVector&>(&World::Evaluate),
             py::arg("e"), py::arg("pl"), py::arg("q"))
        .def("Evaluate", py::overload_cast<ex, ParameterList&, DoubleVector&, std::string, std::string, std::string>(&World::Evaluate),
             py::arg("e"), py::arg("pl"), py::arg("q"), py::arg("fname"), py::arg("comment") = "", py::arg("prefix") = "# ");
    
    // Expose SubUnit classes
    py::class_<SubUnit, std::unique_ptr<SubUnit, py::nodelete>>(m, "SubUnit");
    
    py::class_<GaussianPolymer, SubUnit, std::unique_ptr<GaussianPolymer, py::nodelete>>(m, "GaussianPolymer")
        .def(py::init<>());
    
    py::class_<GaussianLoop, SubUnit, std::unique_ptr<GaussianLoop, py::nodelete>>(m, "GaussianLoop")
        .def(py::init<>());
    
    py::class_<ThinCircle, SubUnit, std::unique_ptr<ThinCircle, py::nodelete>>(m, "ThinCircle")
        .def(py::init<>());
    
    py::class_<ThinRod, SubUnit, std::unique_ptr<ThinRod, py::nodelete>>(m, "ThinRod")
        .def(py::init<>());
    
    py::class_<ThinDisk, SubUnit, std::unique_ptr<ThinDisk, py::nodelete>>(m, "ThinDisk")
        .def(py::init<>());
    
    py::class_<ThinSphericalShell, SubUnit, std::unique_ptr<ThinSphericalShell, py::nodelete>>(m, "ThinSphericalShell")
        .def(py::init<>());
    
    py::class_<SolidSphere, SubUnit, std::unique_ptr<SolidSphere, py::nodelete>>(m, "SolidSphere")
        .def(py::init<>());
    
    py::class_<SolidSphericalShell, SubUnit, std::unique_ptr<SolidSphericalShell, py::nodelete>>(m, "SolidSphericalShell")
        .def(py::init<>());
    
    py::class_<SolidCylinder, SubUnit, std::unique_ptr<SolidCylinder, py::nodelete>>(m, "SolidCylinder")
        .def(py::init<>());
    
    py::class_<Point, SubUnit, std::unique_ptr<Point, py::nodelete>>(m, "Point")
        .def(py::init<>());
    
    // Expose utility functions
    m.def("to_string_python", &to_string_python, "Convert a GiNaC expression to a Python string");
    m.def("to_string_latex", &to_string_latex, "Convert a GiNaC expression to a LaTeX string");
    m.def("to_string_cform", &to_string_cform, "Convert a GiNaC expression to a C/C++ string");
}
