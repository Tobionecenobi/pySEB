#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "SEB.hpp"
#include "World.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "Symbolic.hpp"
#include "SymbolInterface.hpp"
#include "bindingsSymbolic.hpp"
#include "WorldIO/WorldParser.hpp"

namespace py = pybind11;

// Forward declaration of the function that registers common types
void register_common_types(py::module& m);

namespace {

// Transfer ownership of numeric storage without constructing Python floats.
py::array_t<double> numeric_array(DoubleVector values)
{
    auto* storage = new DoubleVector(std::move(values));
    py::capsule owner(storage, [](void* pointer) {
        delete static_cast<DoubleVector*>(pointer);
    });
    // Explicit strides also work with older pybind11 releases on NumPy 2,
    // where deriving strides from the dtype descriptor can yield zero.
    return py::array_t<double>({static_cast<py::ssize_t>(storage->size())},
                               {static_cast<py::ssize_t>(sizeof(double))},
                               storage->data(), owner);
}

} // namespace

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

    py::enum_<IntegrationMethod>(m, "IntegrationMethod")
        .value("QAG", IntegrationMethod::QAG)
        .value("CQUAD", IntegrationMethod::CQUAD);

    py::enum_<LengthUnit>(m, "LengthUnit")
        .value("Meter", LengthUnit::Meter).value("Millimeter", LengthUnit::Millimeter)
        .value("Micrometer", LengthUnit::Micrometer).value("Nanometer", LengthUnit::Nanometer)
        .value("Angstrom", LengthUnit::Angstrom);
    py::enum_<USDLayoutMode>(m, "USDLayoutMode")
        .value("Random", USDLayoutMode::Random)
        .value("Readable", USDLayoutMode::Readable);
    py::class_<USDExportOptions>(m, "USDExportOptions")
        .def(py::init<LengthUnit>())
        .def(py::init<double>())
        .def_readwrite("seed", &USDExportOptions::seed)
        .def_readwrite("layout_mode", &USDExportOptions::layoutMode)
        .def_readwrite("orientation_trials", &USDExportOptions::orientationTrials)
        .def_readwrite("relaxation_sweeps", &USDExportOptions::relaxationSweeps)
        .def_readwrite("minimum_clearance", &USDExportOptions::minimumClearance)
        .def_readwrite("curve_samples", &USDExportOptions::curveSamples)
        .def_readwrite("surface_samples", &USDExportOptions::surfaceSamples)
        .def_readwrite("reference_markers", &USDExportOptions::referenceMarkers)
        .def_readwrite("zero_radius_marker_size", &USDExportOptions::zeroRadiusMarkerSize)
        .def_readwrite("debye_envelope", &USDExportOptions::debyeEnvelope)
        .def_readwrite("debye_envelope_resolution", &USDExportOptions::debyeEnvelopeResolution)
        .def_readwrite("debye_envelope_padding", &USDExportOptions::debyeEnvelopePadding)
        .def_readwrite("debye_envelope_opacity", &USDExportOptions::debyeEnvelopeOpacity)
        .def_readwrite("meters_per_unit", &USDExportOptions::metersPerUnit)
        .def_readwrite("custom_unit", &USDExportOptions::customUnit)
        .def_readwrite("color_overrides", &USDExportOptions::colorOverrides)
        .def_readwrite("opacity_overrides", &USDExportOptions::opacityOverrides);

    py::class_<IntegrationOptions>(m, "SubunitIntegrationOptions")
        .def_readonly("method", &IntegrationOptions::method)
        .def_readonly("absolute_tolerance", &IntegrationOptions::absoluteTolerance)
        .def_readonly("relative_tolerance", &IntegrationOptions::relativeTolerance)
        .def_readonly("workspace_size", &IntegrationOptions::workspaceSize)
        .def_property_readonly("qag_rule", [](const IntegrationOptions& options) {
            static const int points[] = {0, 15, 21, 31, 41, 51, 61};
            return options.qagRule >= 1 && options.qagRule <= 6
                ? points[options.qagRule] : options.qagRule;
        });

    py::class_<pyseb::IntegralDefinition::Dimension>(m, "SubunitIntegralDimension")
        .def_readonly("variable", &pyseb::IntegralDefinition::Dimension::variable)
        .def_property_readonly("lower", [](const pyseb::IntegralDefinition::Dimension& dimension) {
            return dimension.lower.source();
        })
        .def_property_readonly("upper", [](const pyseb::IntegralDefinition::Dimension& dimension) {
            return dimension.upper.source();
        });

    py::class_<pyseb::IntegralDefinition>(m, "SubunitIntegralDefinition")
        .def_readonly("name", &pyseb::IntegralDefinition::name)
        .def_readonly("variable", &pyseb::IntegralDefinition::variable)
        .def_property_readonly("lower", [](const pyseb::IntegralDefinition& integral) {
            return integral.lower.source();
        })
        .def_property_readonly("upper", [](const pyseb::IntegralDefinition& integral) {
            return integral.upper.source();
        })
        .def_property_readonly("integrand", [](const pyseb::IntegralDefinition& integral) {
            return integral.integrand.source();
        })
        .def_readonly("dimensions", &pyseb::IntegralDefinition::dimensions)
        .def_readonly("integration", &pyseb::IntegralDefinition::integration);

    py::class_<pyseb::ParameterDefinition>(m, "SubunitParameterDefinition")
        .def_readonly("name", &pyseb::ParameterDefinition::name)
        .def_readonly("unit", &pyseb::ParameterDefinition::unit)
        .def_readonly("description", &pyseb::ParameterDefinition::description);

    py::class_<pyseb::SubunitMetadata>(m, "SubunitMetadata")
        .def_readonly("title", &pyseb::SubunitMetadata::title)
        .def_readonly("description", &pyseb::SubunitMetadata::description)
        .def_readonly("authors", &pyseb::SubunitMetadata::authors)
        .def_readonly("citations", &pyseb::SubunitMetadata::citations)
        .def_readonly("license", &pyseb::SubunitMetadata::license);

    py::enum_<pyseb::VisualizationGeometryKind>(m, "VisualizationGeometryKind")
        .value("Curve", pyseb::VisualizationGeometryKind::Curve)
        .value("Surface", pyseb::VisualizationGeometryKind::Surface)
        .value("RandomWalk", pyseb::VisualizationGeometryKind::RandomWalk);
    py::class_<pyseb::VisualizationGeometry>(m, "VisualizationGeometry")
        .def_readonly("name", &pyseb::VisualizationGeometry::name)
        .def_readonly("kind", &pyseb::VisualizationGeometry::kind)
        .def_readonly("samples", &pyseb::VisualizationGeometry::samples)
        .def_readonly("distribution", &pyseb::VisualizationGeometry::distribution)
        .def_readonly("closure", &pyseb::VisualizationGeometry::closure);
    py::class_<pyseb::VisualizationReferenceVariant>(m, "VisualizationReferenceVariant")
        .def_readonly("geometry", &pyseb::VisualizationReferenceVariant::geometry)
        .def_readonly("sampling", &pyseb::VisualizationReferenceVariant::sampling);
    py::class_<pyseb::VisualizationReference>(m, "VisualizationReference")
        .def_readonly("name", &pyseb::VisualizationReference::name)
        .def_readonly("kind", &pyseb::VisualizationReference::kind)
        .def_readonly("geometry", &pyseb::VisualizationReference::geometry)
        .def_readonly("sampling", &pyseb::VisualizationReference::sampling)
        .def_readonly("variants", &pyseb::VisualizationReference::variants);
    py::class_<pyseb::VisualizationDefinition>(m, "VisualizationDefinition")
        .def_readonly("present", &pyseb::VisualizationDefinition::present)
        .def_readonly("geometry", &pyseb::VisualizationDefinition::geometry)
        .def_readonly("references", &pyseb::VisualizationDefinition::references)
        .def_readonly("opacity", &pyseb::VisualizationDefinition::opacity)
        .def_readonly("double_sided", &pyseb::VisualizationDefinition::doubleSided)
        .def_readonly("curve_width", &pyseb::VisualizationDefinition::curveWidth);

    py::class_<pyseb::SubunitDefinition>(m, "SubunitDefinition")
        .def_readonly("schema_version", &pyseb::SubunitDefinition::schemaVersion)
        .def_readonly("id", &pyseb::SubunitDefinition::id)
        .def_readonly("api_name", &pyseb::SubunitDefinition::apiName)
        .def_readonly("model_version", &pyseb::SubunitDefinition::modelVersion)
        .def_readonly("source", &pyseb::SubunitDefinition::source)
        .def_readonly("metadata", &pyseb::SubunitDefinition::metadata)
        .def_readonly("invisible", &pyseb::SubunitDefinition::invisible)
        .def_readonly("parameters", &pyseb::SubunitDefinition::parameters)
        .def_readonly("integration", &pyseb::SubunitDefinition::integration)
        .def_readonly("integrals", &pyseb::SubunitDefinition::integrals)
        .def_readonly("specific_references", &pyseb::SubunitDefinition::specificReferences)
        .def_readonly("distributed_references", &pyseb::SubunitDefinition::distributedReferences)
        .def_readonly("visualization", &pyseb::SubunitDefinition::visualization)
        .def_property_readonly("validation_case_count", [](const pyseb::SubunitDefinition& definition) {
            return definition.validationCases.size();
        });

    py::class_<pyseb::SubunitModelInfo>(m, "SubunitModelInfo")
        .def_readonly("id", &pyseb::SubunitModelInfo::id)
        .def_readonly("api_name", &pyseb::SubunitModelInfo::apiName)
        .def_readonly("model_version", &pyseb::SubunitModelInfo::modelVersion)
        .def_readonly("title", &pyseb::SubunitModelInfo::title)
        .def_readonly("source", &pyseb::SubunitModelInfo::source)
        .def_readonly("bundled", &pyseb::SubunitModelInfo::bundled);

    py::class_<pyseb::SubunitValidationFailure>(m, "SubunitValidationFailure")
        .def_readonly("case_name", &pyseb::SubunitValidationFailure::caseName)
        .def_readonly("message", &pyseb::SubunitValidationFailure::message);

    py::class_<pyseb::SubunitValidationReport>(m, "SubunitValidationReport")
        .def_readonly("model_id", &pyseb::SubunitValidationReport::modelId)
        .def_readonly("case_count", &pyseb::SubunitValidationReport::caseCount)
        .def_readonly("failures", &pyseb::SubunitValidationReport::failures)
        .def_readonly("warnings", &pyseb::SubunitValidationReport::warnings)
        .def_property_readonly("ok", &pyseb::SubunitValidationReport::ok);

    m.def("load_subunit_definition", &pyseb::LoadSubunitDefinitionFile, py::arg("path"));
    m.def("validate_subunit_file", &pyseb::ValidateSubunitFile, py::arg("path"));
    m.def("load_world", &pyseb::LoadWorld,
          py::arg("path"), py::arg("model_files") = std::vector<std::string>());
    m.def("world_from_yaml", &pyseb::WorldFromYaml,
          py::arg("yaml"), py::arg("model_files") = std::vector<std::string>());
    m.def("save_world", &pyseb::SaveWorld, py::arg("world"), py::arg("path"));
    m.def("world_to_yaml", &pyseb::WorldToYaml, py::arg("world"));

    py::class_<
        SymbolicSubunit,
        SubUnit,
        std::unique_ptr<SymbolicSubunit, py::nodelete>
    >(m, "SymbolicSubunit")
        .def(py::init<>())
        .def("addReferencePoint", &SymbolicSubunit::setReferencePointName,
             py::arg("reference"));

    py::class_<
        FileDefinedSubunit,
        SubUnit,
        std::unique_ptr<FileDefinedSubunit, py::nodelete>
    >(m, "FileDefinedSubunit")
        .def(py::init<const pyseb::SubunitDefinition&>());

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

    py::class_<CromerMannCoefficients>(m, "CromerMannCoefficients")
        .def(py::init<>())
        .def(py::init<
                 const std::array<double, 4>&,
                 const std::array<double, 4>&,
                 double>(),
             py::arg("a"), py::arg("b"), py::arg("c"))
        .def_readwrite("a", &CromerMannCoefficients::a)
        .def_readwrite("b", &CromerMannCoefficients::b)
        .def_readwrite("c", &CromerMannCoefficients::c)
        .def("evaluate", &CromerMannCoefficients::evaluate, py::arg("q"))
        .def("value_at_zero", &CromerMannCoefficients::valueAtZero)
        .def("second_moment", &CromerMannCoefficients::secondMoment);

    py::class_<AtomicSAXSScatterer>(m, "AtomicSAXSScatterer")
        .def(py::init<>())
        .def_readwrite("center", &AtomicSAXSScatterer::center)
        .def_readwrite("vacuum", &AtomicSAXSScatterer::vacuum)
        .def_readwrite("excluded_volume", &AtomicSAXSScatterer::excludedVolume)
        .def_readwrite("occupancy", &AtomicSAXSScatterer::occupancy)
        .def_readwrite("temperature_factor", &AtomicSAXSScatterer::temperatureFactor)
        .def_readwrite(
            "solvent_accessibility",
            &AtomicSAXSScatterer::solventAccessibility)
        .def_readwrite("implicit_hydrogens", &AtomicSAXSScatterer::implicitHydrogens)
        .def_readwrite("implicit_hydrogen_distance", &AtomicSAXSScatterer::implicitHydrogenDistance)
        .def_readwrite("hydration_volume", &AtomicSAXSScatterer::hydrationVolume)
        .def_readwrite("hydration_uses_solvent_density", &AtomicSAXSScatterer::hydrationUsesSolventDensity)
        .def_readwrite("element", &AtomicSAXSScatterer::element);

    py::enum_<AtomicSAXSHydrationModel>(m, "AtomicSAXSHydrationModel")
        .value("AtomCentered", AtomicSAXSHydrationModel::AtomCentered)
        .value("SpatialShell", AtomicSAXSHydrationModel::SpatialShell);

    py::enum_<AtomicSAXSHydrationWaterModel>(m, "AtomicSAXSHydrationWaterModel")
        .value("NeutralWater", AtomicSAXSHydrationWaterModel::NeutralWater)
        .value("BondedWater", AtomicSAXSHydrationWaterModel::BondedWater)
        .value("FoXSReciprocal", AtomicSAXSHydrationWaterModel::FoXSReciprocal)
        .value("FoXSFast", AtomicSAXSHydrationWaterModel::FoXSFast);

    py::class_<AtomicSAXSOptions>(m, "AtomicSAXSOptions")
        .def(py::init<>())
        .def_readwrite(
            "solvent_electron_density",
            &AtomicSAXSOptions::solventElectronDensity)
        .def_readwrite(
            "excluded_volume_scale",
            &AtomicSAXSOptions::excludedVolumeScale)
        .def_readwrite(
            "include_excluded_solvent",
            &AtomicSAXSOptions::includeExcludedSolvent)
        .def_readwrite(
            "use_temperature_factors",
            &AtomicSAXSOptions::useTemperatureFactors)
        .def_readwrite("include_hydration", &AtomicSAXSOptions::includeHydration)
        .def_readwrite("hydration_model", &AtomicSAXSOptions::hydrationModel)
        .def_readwrite("hydration_water_model", &AtomicSAXSOptions::hydrationWaterModel)
        .def_readwrite(
            "hydration_shell_contrast",
            &AtomicSAXSOptions::hydrationShellContrast)
        .def_readwrite(
            "hydration_volume_per_exposed_atom",
            &AtomicSAXSOptions::hydrationVolumePerExposedAtom);

    py::class_<AtomRecord>(m, "AtomRecord")
        .def(py::init<>())
        .def_readwrite("record_type", &AtomRecord::recordType)
        .def_readwrite("serial", &AtomRecord::serial)
        .def_readwrite("atom_name", &AtomRecord::atomName)
        .def_readwrite("alternate_location", &AtomRecord::alternateLocation)
        .def_readwrite("residue_name", &AtomRecord::residueName)
        .def_readwrite("chain_id", &AtomRecord::chainId)
        .def_readwrite("residue_number", &AtomRecord::residueNumber)
        .def_readwrite("insertion_code", &AtomRecord::insertionCode)
        .def_readwrite("coordinate", &AtomRecord::coordinate)
        .def_readwrite("occupancy", &AtomRecord::occupancy)
        .def_readwrite("temperature_factor", &AtomRecord::temperatureFactor)
        .def_readwrite("element", &AtomRecord::element)
        .def_readwrite("charge", &AtomRecord::charge)
        .def_readwrite("model_number", &AtomRecord::modelNumber);

    py::enum_<AlternateLocationPolicy>(m, "AlternateLocationPolicy")
        .value("Primary", AlternateLocationPolicy::Primary)
        .value("All", AlternateLocationPolicy::All);

    py::class_<StructureParseOptions>(m, "StructureParseOptions")
        .def(py::init<>())
        .def_readwrite("model_number", &StructureParseOptions::modelNumber)
        .def_readwrite("include_hetatm", &StructureParseOptions::includeHetatm)
        .def_readwrite("include_water", &StructureParseOptions::includeWater)
        .def_readwrite("include_hydrogen", &StructureParseOptions::includeHydrogen)
        .def_readwrite(
            "alternate_locations",
            &StructureParseOptions::alternateLocations
        );

    py::class_<StructureParser>(m, "StructureParser");

    py::class_<PDBParser, StructureParser>(m, "PDBParser")
        .def(py::init<>())
        .def("parse_file", &PDBParser::ParseFile,
             py::arg("filename"),
             py::arg("options") = StructureParseOptions())
        .def(
            "parse_string",
            [](const PDBParser& parser,
               const std::string& contents,
               const StructureParseOptions& options) {
                std::istringstream input(contents);
                return parser.Parse(input, options, "<string>");
            },
            py::arg("contents"),
            py::arg("options") = StructureParseOptions()
        );

    py::class_<AtomScattererParameters>(m, "AtomScattererParameters")
        .def(py::init<double, double>(),
             py::arg("radius") = 0.0, py::arg("beta") = 0.0)
        .def_readwrite("radius", &AtomScattererParameters::radius)
        .def_readwrite("beta", &AtomScattererParameters::beta);

    py::class_<AtomParameterProfile>(m, "AtomParameterProfile")
        .def(py::init<>())
        .def(
            "set_element",
            py::overload_cast<const std::string&, double, double>(
                &AtomParameterProfile::SetElement
            ),
            py::arg("element"), py::arg("radius"), py::arg("beta")
        )
        .def(
            "set_atom",
            py::overload_cast<
                const std::string&,
                const std::string&,
                double,
                double
            >(&AtomParameterProfile::SetAtom),
            py::arg("residue_name"), py::arg("atom_name"),
            py::arg("radius"), py::arg("beta")
        )
        .def("has_parameters", &AtomParameterProfile::HasParameters,
             py::arg("atom"))
        .def("resolve", &AtomParameterProfile::Resolve,
             py::arg("atom"));

    py::enum_<AtomicSAXSVolumeModel>(m, "AtomicSAXSVolumeModel")
        .value("BondiElement", AtomicSAXSVolumeModel::BondiElement)
        .value("Fraser1978", AtomicSAXSVolumeModel::Fraser1978)
        .value("Pontius1996", AtomicSAXSVolumeModel::Pontius1996)
        .export_values();

    py::class_<AtomicSAXSParameters>(m, "AtomicSAXSParameters")
        .def(py::init<
                 const CromerMannCoefficients&,
                 double,
                 int,
                 double>(),
             py::arg("vacuum") = CromerMannCoefficients(),
             py::arg("excluded_volume") = 0.0,
             py::arg("implicit_hydrogens") = 0,
             py::arg("sasa_radius") = 0.0)
        .def_readwrite("vacuum", &AtomicSAXSParameters::vacuum)
        .def_readwrite("excluded_volume", &AtomicSAXSParameters::excludedVolume)
        .def_readwrite("implicit_hydrogens", &AtomicSAXSParameters::implicitHydrogens)
        .def_readwrite("atomic_excluded_volume", &AtomicSAXSParameters::atomicExcludedVolume)
        .def_readwrite("implicit_hydrogen_distance", &AtomicSAXSParameters::implicitHydrogenDistance)
        .def_readwrite("sasa_radius", &AtomicSAXSParameters::sasaRadius);

    py::class_<AtomicSAXSParameterProfile>(m, "AtomicSAXSParameterProfile")
        .def(py::init<>())
        .def("set_element", &AtomicSAXSParameterProfile::SetElement,
             py::arg("element"), py::arg("parameters"))
        .def("set_atom", &AtomicSAXSParameterProfile::SetAtom,
             py::arg("residue_name"), py::arg("atom_name"),
             py::arg("parameters"))
        .def("has_parameters", &AtomicSAXSParameterProfile::HasParameters,
             py::arg("atom"))
        .def("resolve", &AtomicSAXSParameterProfile::Resolve,
             py::arg("atom"))
        .def_static(
            "element_defaults",
            &AtomicSAXSParameterProfile::ElementDefaults)
        .def_static(
            "protein_defaults",
            &AtomicSAXSParameterProfile::ProteinDefaults,
            py::arg("volume_model") = AtomicSAXSVolumeModel::Fraser1978);

    py::class_<AtomCloudBuildOptions>(m, "AtomCloudBuildOptions")
        .def(py::init<>())
        .def_readwrite(
            "scale_beta_by_occupancy",
            &AtomCloudBuildOptions::scaleBetaByOccupancy
        )
        .def_readwrite(
            "reference_atom_serials",
            &AtomCloudBuildOptions::referenceAtomSerials
        );

    py::class_<AtomCloudBuilder>(m, "AtomCloudBuilder")
        .def_static(
            "build",
            [](const std::vector<AtomRecord>& atoms,
               const AtomParameterProfile& profile,
               const AtomCloudBuildOptions& options) {
                return AtomCloudBuilder::Build(
                    atoms,
                    profile,
                    options
                ).release();
            },
            py::arg("atoms"),
            py::arg("profile"),
            py::arg("options") = AtomCloudBuildOptions(),
            py::return_value_policy::reference
        );

    py::class_<StructureCloudLoader>(m, "StructureCloudLoader")
        .def_static(
            "load_pdb",
            [](const std::string& filename,
               const AtomParameterProfile& profile,
               const StructureParseOptions& parseOptions,
               const AtomCloudBuildOptions& buildOptions) {
                return StructureCloudLoader::LoadPDB(
                    filename,
                    profile,
                    parseOptions,
                    buildOptions
                ).release();
            },
            py::arg("filename"),
            py::arg("profile"),
            py::arg("parse_options") = StructureParseOptions(),
            py::arg("build_options") = AtomCloudBuildOptions(),
             py::return_value_policy::reference
        );

    py::enum_<AtomicSAXSHydrationShellConvention>(m, "AtomicSAXSHydrationShellConvention")
        .value("VolumeDensity", AtomicSAXSHydrationShellConvention::VolumeDensity)
        .value("Pepsi3", AtomicSAXSHydrationShellConvention::Pepsi3);

    py::class_<AtomicSAXSCloudBuildOptions>(m, "AtomicSAXSCloudBuildOptions")
        .def(py::init<>())
        .def_readwrite("hydration_shell_convention", &AtomicSAXSCloudBuildOptions::hydrationShellConvention)
        .def_readwrite("sasa_use_profile_radii", &AtomicSAXSCloudBuildOptions::sasaUseProfileRadii)
        .def_readwrite("hydration_shell_inner_radius", &AtomicSAXSCloudBuildOptions::hydrationShellInnerRadius)
        .def_readwrite("hydration_shell_thickness", &AtomicSAXSCloudBuildOptions::hydrationShellThickness)
        .def_readwrite("hydration_shell_spacing", &AtomicSAXSCloudBuildOptions::hydrationShellSpacing)
        .def_readwrite(
            "scale_by_occupancy",
            &AtomicSAXSCloudBuildOptions::scaleByOccupancy)
        .def_readwrite(
            "use_implicit_hydrogens",
            &AtomicSAXSCloudBuildOptions::useImplicitHydrogens)
        .def_readwrite(
            "calculate_solvent_accessibility",
            &AtomicSAXSCloudBuildOptions::calculateSolventAccessibility)
        .def_readwrite(
            "sasa_probe_radius",
            &AtomicSAXSCloudBuildOptions::sasaProbeRadius)
        .def_readwrite(
            "sasa_point_count",
            &AtomicSAXSCloudBuildOptions::sasaPointCount)
        .def_readwrite(
            "reference_atom_serials",
            &AtomicSAXSCloudBuildOptions::referenceAtomSerials);

    py::class_<AtomicSAXSCloudBuilder>(m, "AtomicSAXSCloudBuilder")
        .def_static(
            "build",
            [](const std::vector<AtomRecord>& atoms,
               const AtomicSAXSParameterProfile& profile,
               const AtomicSAXSOptions& saxsOptions,
               const AtomicSAXSCloudBuildOptions& buildOptions) {
                return AtomicSAXSCloudBuilder::Build(
                    atoms, profile, saxsOptions, buildOptions).release();
            },
            py::arg("atoms"),
            py::arg("profile"),
            py::arg("saxs_options") = AtomicSAXSOptions(),
            py::arg("build_options") = AtomicSAXSCloudBuildOptions(),
            py::return_value_policy::reference
        );

    py::class_<AtomicSAXSStructureLoader>(m, "AtomicSAXSStructureLoader")
        .def_static(
            "load_pdb",
            [](const std::string& filename,
               const AtomicSAXSParameterProfile& profile,
               const StructureParseOptions& parseOptions,
               const AtomicSAXSOptions& saxsOptions,
               const AtomicSAXSCloudBuildOptions& buildOptions) {
                return AtomicSAXSStructureLoader::LoadPDB(
                    filename,
                    profile,
                    parseOptions,
                    saxsOptions,
                    buildOptions
                ).release();
            },
            py::arg("filename"),
            py::arg("profile"),
            py::arg("parse_options") = StructureParseOptions(),
            py::arg("saxs_options") = AtomicSAXSOptions(),
            py::arg("build_options") = AtomicSAXSCloudBuildOptions(),
            py::return_value_policy::reference
        );

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

    py::class_<
        AtomicSAXSCloud,
        NumericalSubunit,
        std::unique_ptr<AtomicSAXSCloud, py::nodelete>
    >(m, "AtomicSAXSCloud")
        .def(py::init<>())
        .def(py::init<
                 const std::vector<AtomicSAXSScatterer>&,
                 const AtomicSAXSOptions&>(),
             py::arg("scatterers"),
             py::arg("options") = AtomicSAXSOptions())
        .def("add_scatterer", &AtomicSAXSCloud::addScatterer,
             py::arg("scatterer"))
        .def("add_reference_point",
             py::overload_cast<const refPoint&, const CartesianPoint3D&>(
                 &AtomicSAXSCloud::addReferencePoint),
             py::arg("name"), py::arg("coordinate"))
        .def("add_reference_point",
             py::overload_cast<const refPoint&, double, double, double>(
                 &AtomicSAXSCloud::addReferencePoint),
             py::arg("name"), py::arg("x"), py::arg("y"), py::arg("z"))
        .def("get_scatterers", &AtomicSAXSCloud::getScatterers,
             py::return_value_policy::reference_internal)
        .def("get_reference_coordinates",
             &AtomicSAXSCloud::getReferenceCoordinates)
        .def("scatterer_count", &AtomicSAXSCloud::scattererCount)
        .def("get_options", &AtomicSAXSCloud::getOptions)
        .def("set_options", &AtomicSAXSCloud::setOptions, py::arg("options"))
        .def("set_solvent_electron_density",
             &AtomicSAXSCloud::setSolventElectronDensity, py::arg("value"))
        .def("set_excluded_volume_scale",
             &AtomicSAXSCloud::setExcludedVolumeScale, py::arg("value"))
        .def("set_include_excluded_solvent",
             &AtomicSAXSCloud::setIncludeExcludedSolvent, py::arg("value"))
        .def("set_use_temperature_factors",
             &AtomicSAXSCloud::setUseTemperatureFactors, py::arg("value"))
        .def("set_hydration", &AtomicSAXSCloud::setHydration,
             py::arg("enabled"), py::arg("shell_contrast") = 0.0,
             py::arg("volume_per_exposed_atom") = 30.0)
        .def("atomic_amplitude", &AtomicSAXSCloud::atomicAmplitude,
             py::arg("index"), py::arg("q"))
        .def_static("sinc", &AtomicSAXSCloud::sinc)
        .def_static(
            "hydrogen_coefficients",
            &AtomicSAXSCloud::hydrogenCoefficients);

    // Expose World class - basic structure only, symbolic methods are registered in backend-specific files
    py::class_<World> world(m, "World");
    world.def(py::init<std::string>(), py::arg("id") = "World")
        .def("to_yaml", &pyseb::WorldToYaml)
        .def("save", [](const World& self, const std::string& path) {
            pyseb::SaveWorld(self, path);
        }, py::arg("path"))
        .def("register_subunit_file", &World::RegisterSubunitFile, py::arg("path"))
        .def("register_subunit_directory", &World::RegisterSubunitDirectory, py::arg("path"))
        .def("list_subunit_models", &World::ListSubunitModels)
        .def("add_subunit_file", &World::AddFile,
             py::arg("path"), py::arg("name"), py::arg("tag") = "")
        .def("add_file", &World::AddFile,
             py::arg("path"), py::arg("name"), py::arg("tag") = "")
        .def("link_subunit_file", &World::LinkFile,
             py::arg("path"), py::arg("new_reference"), py::arg("old_reference"), py::arg("tag") = "")
        .def("link_file", &World::LinkFile,
             py::arg("path"), py::arg("new_reference"), py::arg("old_reference"), py::arg("tag") = "")
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
        .def("_evaluate_form_factor_numpy", [](World& self, const std::string& name,
                const ParameterList& parameters,
                py::array_t<double, py::array::c_style | py::array::forcecast> q,
                const std::string& reference) -> py::object {
            if (q.ndim() != 1) throw py::value_error("Expected a 1-D q array");
            // A numeric copy preserves the C++ vector API without routing
            // every input/output sample through a boxed Python float.
            DoubleVector grid(static_cast<std::size_t>(q.size()));
            std::copy_n(q.data(), q.size(), grid.begin());
            if (reference.empty()) {
                return numeric_array(self.EvaluateFormFactorUnnormalized(name, parameters, grid));
            }
            auto result = self.EvaluateFormFactorAndAmplitudeUnnormalized(
                name, reference, parameters, grid);
            return py::make_tuple(numeric_array(std::move(result.first)),
                                  numeric_array(std::move(result.second)));
        }, py::arg("name"), py::arg("parameters"), py::arg("q"),
           py::arg("reference") = "")
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
             py::arg("parameters"))
        .def("export_usd", &World::ExportUSD,
             py::arg("structure"), py::arg("path"), py::arg("parameters"), py::arg("options"));

    register_symbolic_world_bindings(world);
}
