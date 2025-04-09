cmake_minimum_required(VERSION 3.18)
project(pyseb)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Find required libraries
find_package(GSL REQUIRED)
find_library(GINAC_LIBRARY NAMES ginac)
find_library(CLN_LIBRARY NAMES cln)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/SEB)

# Add the SEB library
add_library(seb STATIC
    SEB/Abstract_subunit.hpp
    SEB/Constants.hpp
    SEB/Exceptions.hpp
    SEB/SEB.hpp
    SEB/SpecialFunctions.cpp
    SEB/SpecialFunctions.hpp
    SEB/Structure.cpp
    SEB/Structure.hpp
    SEB/Subunit.cpp
    SEB/Subunit.hpp
    SEB/SymbolInterface.cpp
    SEB/SymbolInterface.hpp
    SEB/Types.hpp
    SEB/World.cpp
    SEB/World.hpp
    SEB/Subunits/GaussianLoop.cpp
    SEB/Subunits/GaussianPolymer.cpp
    SEB/Subunits/Point.cpp
    SEB/Subunits/SolidCylinder.cpp
    SEB/Subunits/SolidSphere.cpp
    SEB/Subunits/SolidSphericalShell.cpp
    SEB/Subunits/Subunits.hpp
    SEB/Subunits/Symbolic.cpp
    SEB/Subunits/ThinCircle.cpp
    SEB/Subunits/ThinDisk.cpp
    SEB/Subunits/ThinRod.cpp
    SEB/Subunits/ThinSphericalShell.cpp
)

# Link libraries to SEB
target_link_libraries(seb PRIVATE
    ${GINAC_LIBRARY}
    ${CLN_LIBRARY}
    GSL::gsl
    GSL::gslcblas
)

# Create the Python module
pybind11_add_module(_pyseb src/bindings.cpp)
target_link_libraries(_pyseb PRIVATE seb)

# Install the Python module
install(TARGETS _pyseb DESTINATION .)
