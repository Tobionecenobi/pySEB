#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "SEB.hpp"
#include "World.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"
#include "SymbolicInterface.hpp"

namespace py = pybind11;

// Define constants if not already defined
#ifndef GENERIC
#define GENERIC 0
#endif

#ifndef XVAR
#define XVAR 1
#endif

#ifndef QVAR
#define QVAR 2
#endif

// This function registers all the common types that are used in multiple binding files
void register_common_types(py::module& m) {
    // refPoint/subName/structName are typedefs of std::string, and GraphID is typedef int.
    // Python should pass str/int directly for these aliases.

    // Expose SubUnit class
    py::class_<SubUnit>(m, "SubUnit")
        .def("getType", &SubUnit::getType)
        .def("getName", &SubUnit::getName)
        .def("getTag", &SubUnit::getTag)
        .def("getReferencePoints", &SubUnit::getReferencePoints);

    // Expose constants
    m.attr("GENERIC") = GENERIC;
    m.attr("XVAR") = XVAR;
    m.attr("QVAR") = QVAR;
    m.attr("WORLDMAXDEPTH") = WORLDMAXDEPTH;
}
