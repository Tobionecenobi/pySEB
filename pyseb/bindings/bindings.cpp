#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <algorithm>
#include <memory>
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

    py::enum_<NormalizationMode>(m, "NormalizationMode")
        .value("Normalized", NormalizationMode::Normalized)
        .value("Unnormalized", NormalizationMode::Unnormalized);

    py::class_<
        SymbolicSubunit,
        SubUnit,
        std::unique_ptr<SymbolicSubunit, py::nodelete>
    >(m, "SymbolicSubunit")
        .def(py::init<>())
        .def("addReferencePoint", &SymbolicSubunit::setReferencePointName,
             py::arg("reference"));

    py::class_<
        NumericalSubunit,
        SymbolicSubunit,
        std::unique_ptr<NumericalSubunit, py::nodelete>
    >(m, "NumericalSubunit")
        .def(py::init<NormalizationMode>(),
             py::arg("normalization_mode") = NormalizationMode::Normalized)
        .def("addReferencePoint", &NumericalSubunit::addReferencePoint,
             py::arg("reference"))
        .def("addDistributedReferencePointType",
             &NumericalSubunit::addDistributedReferencePointType,
             py::arg("reference"))
        .def("setTotalBeta", &NumericalSubunit::setTotalBeta,
             py::arg("beta"))
        .def("setTotalBetaProvider", &NumericalSubunit::setTotalBetaProvider,
             py::arg("provider"))
        .def("setFormFactorFunction", &NumericalSubunit::setFormFactorFunction,
             py::arg("function"))
        .def("setFormFactorAmplitudeFunction",
             &NumericalSubunit::setFormFactorAmplitudeFunction,
             py::arg("reference"), py::arg("function"))
        .def("setPhaseFactorFunction",
             &NumericalSubunit::setPhaseFactorFunction,
             py::arg("reference1"), py::arg("reference2"), py::arg("function"))
        .def("setRadiusOfGyration2", &NumericalSubunit::setRadiusOfGyration2,
             py::arg("value"))
        .def("setRadiusOfGyration2Provider",
             &NumericalSubunit::setRadiusOfGyration2Provider,
             py::arg("provider"))
        .def("setSigmaMSDRef2Scat", &NumericalSubunit::setSigmaMSDRef2Scat,
             py::arg("reference"), py::arg("value"))
        .def("setSigmaMSDRef2ScatProvider",
             &NumericalSubunit::setSigmaMSDRef2ScatProvider,
             py::arg("reference"), py::arg("provider"))
        .def("setSigmaMSDRef2Ref", &NumericalSubunit::setSigmaMSDRef2Ref,
             py::arg("reference1"), py::arg("reference2"), py::arg("value"))
        .def("setSigmaMSDRef2RefProvider",
             &NumericalSubunit::setSigmaMSDRef2RefProvider,
             py::arg("reference1"), py::arg("reference2"), py::arg("provider"))
        .def("ValidateNumerically", &NumericalSubunit::ValidateNumerically,
             py::arg("parameters") = ParameterList(),
             py::arg("tolerance") = 1e-8)
        .def("getNormalizationMode", &NumericalSubunit::getNormalizationMode);

    py::class_<CartesianPoint3D>(m, "CartesianPoint3D")
        .def(py::init<double, double, double>(),
             py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0)
        .def_readwrite("x", &CartesianPoint3D::x)
        .def_readwrite("y", &CartesianPoint3D::y)
        .def_readwrite("z", &CartesianPoint3D::z);

    py::class_<SphereScatterer>(m, "SphereScatterer")
        .def(py::init<double, double, double, double, double>(),
             py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("radius"), py::arg("beta"))
        .def_readwrite("center", &SphereScatterer::center)
        .def_readwrite("radius", &SphereScatterer::radius)
        .def_readwrite("beta", &SphereScatterer::beta);

    py::class_<
        DebyeSphereCloud,
        NumericalSubunit,
        std::unique_ptr<DebyeSphereCloud, py::nodelete>
    >(m, "DebyeSphereCloud")
        .def(py::init<>())
        .def(py::init<const std::vector<SphereScatterer>&>(),
             py::arg("scatterers"))
        .def("addScatterer",
             py::overload_cast<const SphereScatterer&>(
                 &DebyeSphereCloud::addScatterer),
             py::arg("scatterer"))
        .def("addScatterer",
             py::overload_cast<double, double, double, double, double>(
                 &DebyeSphereCloud::addScatterer),
             py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("radius"), py::arg("beta"))
        .def("addReferencePoint",
             py::overload_cast<const refPoint&, const CartesianPoint3D&>(
                 &DebyeSphereCloud::addReferencePoint),
             py::arg("name"), py::arg("coordinate"))
        .def("addReferencePoint",
             py::overload_cast<const refPoint&, double, double, double>(
                 &DebyeSphereCloud::addReferencePoint),
             py::arg("name"), py::arg("x"), py::arg("y"), py::arg("z"))
        .def("getScatterers", &DebyeSphereCloud::getScatterers,
             py::return_value_policy::reference_internal)
        .def("getReferenceCoordinates",
             &DebyeSphereCloud::getReferenceCoordinates)
        .def("scattererCount", &DebyeSphereCloud::scattererCount)
        .def_static("sinc", &DebyeSphereCloud::sinc)
        .def_static("sphereAmplitude", &DebyeSphereCloud::sphereAmplitude);

    py::class_<GaussianPolymer, SubUnit, std::unique_ptr<GaussianPolymer, py::nodelete>>(
        m, "GaussianPolymer"
    ).def(py::init<>());
    py::class_<GaussianLoop, SubUnit, std::unique_ptr<GaussianLoop, py::nodelete>>(
        m, "GaussianLoop"
    ).def(py::init<>());
    py::class_<ThinCircle, SubUnit, std::unique_ptr<ThinCircle, py::nodelete>>(
        m, "ThinCircle"
    ).def(py::init<>());
    py::class_<ThinRod, SubUnit, std::unique_ptr<ThinRod, py::nodelete>>(
        m, "ThinRod"
    ).def(py::init<>());
    py::class_<ThinDisk, SubUnit, std::unique_ptr<ThinDisk, py::nodelete>>(
        m, "ThinDisk"
    ).def(py::init<>());
    py::class_<
        ThinSphericalShell,
        SubUnit,
        std::unique_ptr<ThinSphericalShell, py::nodelete>
    >(m, "ThinSphericalShell").def(py::init<>());
    py::class_<SolidSphere, SubUnit, std::unique_ptr<SolidSphere, py::nodelete>>(
        m, "SolidSphere"
    ).def(py::init<>());
    py::class_<
        SolidSphericalShell,
        SubUnit,
        std::unique_ptr<SolidSphericalShell, py::nodelete>
    >(m, "SolidSphericalShell").def(py::init<>());
    py::class_<SolidCylinder, SubUnit, std::unique_ptr<SolidCylinder, py::nodelete>>(
        m, "SolidCylinder"
    ).def(py::init<>());
    py::class_<Point, SubUnit, std::unique_ptr<Point, py::nodelete>>(
        m, "Point"
    ).def(py::init<>());

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
             py::arg("sub"), py::arg("name"), py::arg("tag") = "",
             py::keep_alive<1, 2>())
        .def("Add", py::overload_cast<std::string, subName, std::string>(&World::Add),
             py::arg("subunit_type"), py::arg("name"), py::arg("tag") = "")
        .def("Add", py::overload_cast<GraphID, structName>(&World::Add),
             py::arg("gid"), py::arg("name"))
        .def("Link", py::overload_cast<SubUnit*, refPoint, refPoint, std::string>(&World::Link),
             py::arg("sub"), py::arg("newr"), py::arg("oldr"), py::arg("tag") = "",
             py::keep_alive<1, 2>())
        .def("Link", py::overload_cast<std::string, refPoint, refPoint, std::string>(&World::Link),
             py::arg("subunit_type"), py::arg("newr"), py::arg("oldr"), py::arg("tag") = "")
        .def("Link", py::overload_cast<GraphID, refPoint, refPoint>(&World::Link),
             py::arg("gid"), py::arg("r1"), py::arg("r2"))
        .def("getParams", &World::getParams)
        .def("getParamsq", &World::getParamsq)
        .def("getq", &World::getq, py::arg("value") = 0)
        .def("q", &World::q, py::arg("value") = 0)
        .def("EvaluateFormFactor",
             py::overload_cast<std::string, const ParameterList&, double>(
                 &World::EvaluateFormFactor),
             py::arg("name"), py::arg("parameters"), py::arg("q"))
        .def("EvaluateFormFactor",
             py::overload_cast<std::string, const ParameterList&, const DoubleVector&>(
                 &World::EvaluateFormFactor),
             py::arg("name"), py::arg("parameters"), py::arg("q"))
        .def("EvaluateFormFactorUnnormalized",
             py::overload_cast<std::string, const ParameterList&, double>(
                 &World::EvaluateFormFactorUnnormalized),
             py::arg("name"), py::arg("parameters"), py::arg("q"))
        .def("EvaluateFormFactorUnnormalized",
             py::overload_cast<std::string, const ParameterList&, const DoubleVector&>(
                 &World::EvaluateFormFactorUnnormalized),
             py::arg("name"), py::arg("parameters"), py::arg("q"))
        .def("EvaluateFormFactorAmplitude",
             py::overload_cast<refPoint, const ParameterList&, double>(
                 &World::EvaluateFormFactorAmplitude),
             py::arg("reference"), py::arg("parameters"), py::arg("q"))
        .def("EvaluateFormFactorAmplitude",
             py::overload_cast<refPoint, const ParameterList&, const DoubleVector&>(
                 &World::EvaluateFormFactorAmplitude),
             py::arg("reference"), py::arg("parameters"), py::arg("q"))
        .def("EvaluateFormFactorAmplitudeUnnormalized",
             py::overload_cast<refPoint, const ParameterList&, double>(
                 &World::EvaluateFormFactorAmplitudeUnnormalized),
             py::arg("reference"), py::arg("parameters"), py::arg("q"))
        .def("EvaluateFormFactorAmplitudeUnnormalized",
             py::overload_cast<refPoint, const ParameterList&, const DoubleVector&>(
                 &World::EvaluateFormFactorAmplitudeUnnormalized),
             py::arg("reference"), py::arg("parameters"), py::arg("q"))
        .def("EvaluatePhaseFactor",
             py::overload_cast<refPoint, refPoint, const ParameterList&, double>(
                 &World::EvaluatePhaseFactor),
             py::arg("reference1"), py::arg("reference2"),
             py::arg("parameters"), py::arg("q"))
        .def("EvaluatePhaseFactor",
             py::overload_cast<refPoint, refPoint, const ParameterList&, const DoubleVector&>(
                 &World::EvaluatePhaseFactor),
             py::arg("reference1"), py::arg("reference2"),
             py::arg("parameters"), py::arg("q"))
        .def("EvaluateRadiusOfGyration2", &World::EvaluateRadiusOfGyration2,
             py::arg("name"), py::arg("parameters"))
        .def("EvaluateSMSDRef2Scat", &World::EvaluateSMSDRef2Scat,
             py::arg("reference"), py::arg("parameters"))
        .def("EvaluateSMSDRef2Ref", &World::EvaluateSMSDRef2Ref,
             py::arg("reference1"), py::arg("reference2"),
             py::arg("parameters"));

    register_symbolic_world_bindings(world);
}
