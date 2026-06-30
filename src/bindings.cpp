#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <algorithm>
#include <stdexcept>

#include "SEB.hpp"
#include "World.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "Symbolic.hpp"
#include "SymbolInterface.hpp"
#include "bindingsSymbolic.hpp"

namespace py = pybind11;

// Forward declaration of the function that registers common types
void register_common_types(py::module& m);

PYBIND11_MODULE(_pyseb, m) {
    m.doc() = "Python bindings for the Scattering Equation Builder (SEB) library";

    sebsym::initialize();
    SymbolicFactory::instance();

    m.def("available_backends", []() {
        auto backends = SymbolicFactory::availableBackends();
        if (std::find(backends.begin(), backends.end(), "sympy") == backends.end()) {
            backends.push_back("sympy");
        }
        return backends;
    });
    m.def("get_backend", []() {
        return SymbolicFactory::activeBackendName();
    });
    m.def("set_backend", [](const std::string& name) {
        if (name == "sympy") {
            if (!SymbolicFactory::setBackend("portable")) {
                throw std::runtime_error("portable backend is not available for SymPy adapter");
            }
            return std::string("sympy");
        }
        if (!SymbolicFactory::setBackend(name)) {
            throw std::runtime_error("Unknown symbolic backend: " + name);
        }
        return SymbolicFactory::activeBackendName();
    });

    // Register common types
    register_common_types(m);

    // Expose World class - basic structure only, symbolic methods are registered in backend-specific files
    py::class_<World> world(m, "World");
    world.def(py::init<std::string>(), py::arg("id") = "World")
        .def("Add", [](World& self, const std::string& subunit_type) {
            return self.Add(subunit_type, subunit_type);
        }, py::arg("subunit_type"))
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

    register_symbolic_world_bindings(world);
}
