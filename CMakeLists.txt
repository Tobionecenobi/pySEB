cmake_minimum_required(VERSION 3.18)
project(pyseb VERSION 0.1.0)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add the cmake directory to the module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Option to build C++ code with GiNaC (default OFF for Python usage)
option(USE_GINAC "Use GiNaC for C++ symbolic computations" OFF)

# Option to build with SymPy support
option(USE_SYMPY "Use SymPy for Python symbolic computations" ON)

# Option to use GiNaC implementation
option(USE_GINAC_IMPL "Use GiNaC implementation for symbolic interface" OFF)

if(USE_SYMPY)
    message(STATUS "Building with SymPy support (for Python)")
endif()

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Find required libraries
find_package(GSL REQUIRED)

if(USE_GINAC)
    find_library(GINAC_LIBRARY NAMES ginac)
    find_library(CLN_LIBRARY NAMES cln)
    if(GINAC_LIBRARY AND CLN_LIBRARY)
        message(STATUS "Building with GiNaC support")
    else()
        message(WARNING "GiNaC requested but not found, building without GiNaC")
        set(USE_GINAC OFF)
    endif()
endif()

# Add the reusable symbolic backend library
add_library(seb-symbolic STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB/SymbolicInterface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB/SymbolicPortable.cpp
)

target_include_directories(seb-symbolic PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SEB>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(seb-symbolic PUBLIC
    GSL::gsl
    GSL::gslcblas
)

target_compile_definitions(seb-symbolic PUBLIC SPECIALFUNCTIONS_IMPL)

set_target_properties(seb-symbolic PROPERTIES
    EXPORT_NAME symbolic
    PUBLIC_HEADER
        "${CMAKE_CURRENT_SOURCE_DIR}/SEB/Symbolic.hpp;${CMAKE_CURRENT_SOURCE_DIR}/SEB/SymbolicInterface.hpp;${CMAKE_CURRENT_SOURCE_DIR}/SEB/Expression.hpp"
)

add_library(SEB::symbolic ALIAS seb-symbolic)

# Add the SEB library
add_library(seb STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB/SpecialFunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB/SpecialFunctionsWrapper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB/Subunit.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB/SymbolInterface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB/World.cpp
)

if(USE_GINAC)
    target_sources(seb-symbolic PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/SEB/GiNaCSymbolic.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SEB/RegisterSpecialFunctions.cpp
    )

    # Add GiNaC-specific include directories and libraries
    find_package(GiNaC REQUIRED)
    find_package(CLN REQUIRED)
    target_include_directories(seb-symbolic PUBLIC ${GINAC_INCLUDE_DIRS})
    target_link_libraries(seb-symbolic PUBLIC ${GINAC_LIBRARIES} ${CLN_LIBRARY})

    target_compile_definitions(seb-symbolic PUBLIC USE_GINAC USE_GINAC_IMPL)
else()
    # When USE_GINAC is OFF, we're using SymPy
    message(STATUS "Building with SymPy support only (no GiNaC)")

    target_compile_definitions(seb-symbolic PUBLIC USE_SYMPY_IMPL)
endif()

# Include directories for SEB library
target_include_directories(seb PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB
    ${Python_INCLUDE_DIRS}
    ${pybind11_INCLUDE_DIRS}
)

# Link libraries to SEB
target_link_libraries(seb PUBLIC
    seb-symbolic
)

if(USE_SYMPY)
    target_compile_definitions(seb PUBLIC USE_SYMPY)
endif()

if(USE_GINAC)
    target_link_libraries(seb PUBLIC
        ${GINAC_LIBRARY}
        ${CLN_LIBRARY}
    )
endif()

# Create the Python module
set(BINDING_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bindingsTypes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/bindingsSymbolic.cpp
)

pybind11_add_module(_pyseb MODULE ${BINDING_SOURCES})
target_link_libraries(_pyseb PRIVATE seb)

if(USE_GINAC)
    target_compile_definitions(_pyseb PRIVATE USE_GINAC_BINDINGS)
endif()

# Ensure Python bindings can find headers
target_include_directories(_pyseb PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/SEB
)

# Set output directory for both build and install
set_target_properties(_pyseb PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pyseb
    INSTALL_RPATH_USE_LINK_PATH ON
    BUILD_WITH_INSTALL_RPATH ON
    INSTALL_RPATH "$ORIGIN"
)

# Install the Python module into pyseb package directory
install(TARGETS _pyseb DESTINATION pyseb)
install(TARGETS seb-symbolic
        EXPORT SEBSymbolicTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SEB)

install(EXPORT SEBSymbolicTargets
        FILE SEBSymbolicTargets.cmake
        NAMESPACE SEB::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SEBSymbolic)

set(SEB_SYMBOLIC_WITH_GINAC ${USE_GINAC})
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SEBSymbolicConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SEBSymbolicConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SEBSymbolic
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SEBSymbolicConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SEBSymbolicConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/SEBSymbolicConfigVersion.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGiNaC.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindCLN.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SEBSymbolic)

# Create Python examples directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/python)

# Install Python examples
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/python/
        DESTINATION examples/python
        FILES_MATCHING PATTERN "*.py")
