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

# Build options
option(BUILD_PYTHON "Build Python bindings" OFF)
option(SEB_ENABLE_PORTABLE_BACKEND "Enable the portable C++ symbolic backend" ON)
option(SEB_ENABLE_GINAC_BACKEND "Enable the GiNaC C++ symbolic backend" OFF)

# Temporary compatibility aliases for older build commands.
if(DEFINED USE_GINAC AND USE_GINAC)
    set(SEB_ENABLE_GINAC_BACKEND ON CACHE BOOL "Enable the GiNaC C++ symbolic backend" FORCE)
endif()

if(DEFINED USE_GINAC_IMPL AND USE_GINAC_IMPL)
    set(SEB_ENABLE_GINAC_BACKEND ON CACHE BOOL "Enable the GiNaC C++ symbolic backend" FORCE)
endif()

set(USE_GINAC ${SEB_ENABLE_GINAC_BACKEND} CACHE BOOL "Compatibility alias for SEB_ENABLE_GINAC_BACKEND" FORCE)
set(USE_GINAC_IMPL ${SEB_ENABLE_GINAC_BACKEND} CACHE BOOL "Compatibility alias for SEB_ENABLE_GINAC_BACKEND" FORCE)

if(NOT SEB_ENABLE_PORTABLE_BACKEND)
    message(FATAL_ERROR "The portable symbolic backend is required and must remain enabled.")
endif()

# Find required libraries
find_package(GSL REQUIRED)

if(SEB_ENABLE_GINAC_BACKEND)
    message(STATUS "Building with GiNaC symbolic backend")
    find_package(GiNaC REQUIRED)
    find_package(CLN REQUIRED)
endif()

# Add the reusable symbolic backend library
add_library(seb-symbolic STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic/SymbolicInterface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic/SymbolicPortable.cpp
)

target_include_directories(seb-symbolic PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/seb>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/seb-symbolic>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/seb>
)

target_link_libraries(seb-symbolic PUBLIC
    GSL::gsl
    GSL::gslcblas
)

target_compile_definitions(seb-symbolic PUBLIC SPECIALFUNCTIONS_IMPL)

set_target_properties(seb-symbolic PROPERTIES EXPORT_NAME symbolic)

add_library(SEB::symbolic ALIAS seb-symbolic)

# Add the SEB library
add_library(seb STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SpecialFunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SpecialFunctionsWrapper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/Subunit.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SymbolInterface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/World.cpp
)

if(SEB_ENABLE_GINAC_BACKEND)
    target_sources(seb-symbolic PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic/GiNaCSymbolic.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic/RegisterSpecialFunctions.cpp
    )

    target_include_directories(seb-symbolic PUBLIC ${GINAC_INCLUDE_DIRS})
    target_link_libraries(seb-symbolic PUBLIC ${GINAC_LIBRARIES} ${CLN_LIBRARIES})

    target_compile_definitions(seb-symbolic PUBLIC USE_GINAC USE_GINAC_IMPL)
endif()

# Include directories for SEB library
target_include_directories(seb PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/seb>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/seb>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/seb-symbolic>
)

# Link libraries to SEB
target_link_libraries(seb PUBLIC
    seb-symbolic
)

set_target_properties(seb PROPERTIES
    EXPORT_NAME seb
)

add_library(SEB::seb ALIAS seb)

if(BUILD_PYTHON)
    message(STATUS "Building Python bindings")

    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    find_package(pybind11 CONFIG REQUIRED)

    set(BINDING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/bindings/bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/bindings/bindingsTypes.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/bindings/bindingsSymbolic.cpp
    )

    pybind11_add_module(_pyseb MODULE ${BINDING_SOURCES})
    target_link_libraries(_pyseb PRIVATE SEB::seb)

    if(SEB_ENABLE_GINAC_BACKEND)
        target_compile_definitions(_pyseb PRIVATE USE_GINAC_BINDINGS)
    endif()

    target_include_directories(_pyseb PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/seb
        ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic
        ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/bindings
    )

    set_target_properties(_pyseb PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pyseb
        INSTALL_RPATH_USE_LINK_PATH ON
        BUILD_WITH_INSTALL_RPATH ON
        INSTALL_RPATH "$ORIGIN"
    )

    install(TARGETS _pyseb DESTINATION pyseb)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/examples/python)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/examples/python/
            DESTINATION examples/python
            FILES_MATCHING PATTERN "*.py")
endif()

install(TARGETS seb seb-symbolic
        EXPORT SEBSymbolicTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/seb/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/seb
        FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "build" EXCLUDE
            PATTERN "__pycache__" EXCLUDE)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/seb-symbolic
        FILES_MATCHING
            PATTERN "*.hpp"
            PATTERN "build" EXCLUDE
            PATTERN "__pycache__" EXCLUDE)

install(EXPORT SEBSymbolicTargets
        FILE SEBSymbolicTargets.cmake
        NAMESPACE SEB::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SEBSymbolic)

set(SEB_SYMBOLIC_WITH_GINAC ${SEB_ENABLE_GINAC_BACKEND})
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
