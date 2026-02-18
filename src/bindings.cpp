#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "SEB.hpp"
#include "World.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "SymbolInterface.hpp"

namespace py = pybind11;

// Forward declaration of the function that registers common types
void register_common_types(py::module& m);

// Include the appropriate bindings based on the build configuration
#ifdef USE_GINAC_BINDINGS
#include "bindingsGiNaC.hpp"
#include "GiNaCSymbolic.hpp"
#endif

#ifdef USE_SYMPY
#include "bindingsSymPy.hpp"
#else
#include "bindingsFallback.hpp"
#endif

PYBIND11_MODULE(_pyseb, m) {
    m.doc() = "Python bindings for the Scattering Equation Builder (SEB) library";

#ifdef USE_GINAC_BINDINGS
    static GiNaCFactory ginac_factory;
    SymbolicFactory::setInstance(&ginac_factory);
#endif

    // Register common types
    register_common_types(m);

    // Expose World class - basic structure only, symbolic methods are registered in backend-specific files
    py::class_<World> world(m, "World");
    world.def(py::init<std::string>(), py::arg("id") = "World")
        .def("Add", [](World& self, const std::string& subunit_type, const std::string& tag) {
            return self.Add(subunit_type, subunit_type, tag);
        }, py::arg("subunit_type"), py::arg("tag") = "")
        .def("Add", [](World& self, GraphID gid, const std::string& name) {
            return self.Add(gid, name);
        }, py::arg("gid"), py::arg("name"))
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
        .def("getParams", &World::getParams)
        .def("getParamsq", &World::getParamsq)
        .def("getq", &World::getq, py::arg("value") = 0)
        .def("q", &World::q, py::arg("value") = 0);

    // Register backend-specific bindings
#ifdef USE_GINAC_BINDINGS
    register_ginac_bindings(world);
#endif

#ifdef USE_SYMPY
    register_sympy_bindings(world);
#else
    register_fallback_bindings(world);
#endif
}
