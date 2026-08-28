cmake_minimum_required(VERSION 3.18)
project(pyseb VERSION 0.2.2)

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
option(BUILD_DOC_TUTORIALS "Build C++ programs used by the documentation" OFF)
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

# Bundle every immediate model file into the core library so adding a catalogue
# model requires only adding its YAML definition. CONFIGURE_DEPENDS makes CMake
# regenerate the catalogue when files are added or removed.
file(GLOB PYSEB_BUNDLED_MODEL_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/models/*.pyseb.yaml"
)
list(SORT PYSEB_BUNDLED_MODEL_FILES)
if(NOT PYSEB_BUNDLED_MODEL_FILES)
    message(FATAL_ERROR "No bundled models found in ${CMAKE_CURRENT_SOURCE_DIR}/models")
endif()

set(PYSEB_BUNDLED_MODEL_ENTRIES "")
foreach(PYSEB_MODEL_FILE IN LISTS PYSEB_BUNDLED_MODEL_FILES)
    get_filename_component(PYSEB_MODEL_FILENAME "${PYSEB_MODEL_FILE}" NAME)
    if(NOT PYSEB_MODEL_FILENAME MATCHES "^[A-Za-z0-9_.-]+$")
        message(FATAL_ERROR "Bundled model filename contains unsupported characters: ${PYSEB_MODEL_FILENAME}")
    endif()
    file(READ "${PYSEB_MODEL_FILE}" PYSEB_MODEL_YAML)
    string(MD5 PYSEB_MODEL_HASH "${PYSEB_MODEL_FILENAME};${PYSEB_MODEL_YAML}")
    string(SUBSTRING "${PYSEB_MODEL_HASH}" 0 10 PYSEB_MODEL_HASH_SHORT)
    set(PYSEB_MODEL_DELIMITER "PY${PYSEB_MODEL_HASH_SHORT}")
    string(FIND "${PYSEB_MODEL_YAML}" ")${PYSEB_MODEL_DELIMITER}\"" PYSEB_DELIMITER_COLLISION)
    if(NOT PYSEB_DELIMITER_COLLISION EQUAL -1)
        message(FATAL_ERROR "Unable to embed ${PYSEB_MODEL_FILENAME}: raw-string delimiter collision")
    endif()
    string(APPEND PYSEB_BUNDLED_MODEL_ENTRIES
        "    {\"${PYSEB_MODEL_FILENAME}\", R\"${PYSEB_MODEL_DELIMITER}(${PYSEB_MODEL_YAML})${PYSEB_MODEL_DELIMITER}\"},\n")
endforeach()
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SubunitIO/BundledModels.cpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/BundledModels.cpp
    @ONLY
)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/NumericalIntegration.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SpecialFunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SpecialFunctionsWrapper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/Subunit.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SubunitIO/ExpressionParser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SubunitIO/SubunitDefinition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SubunitIO/SubunitRegistry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SubunitIO/SubunitValidator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/Subunits/FileDefined.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/generated/BundledModels.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/SymbolInterface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/World.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/USDExport.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/StructureIO/PDBParser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/StructureIO/AtomParameterProfile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/seb/StructureIO/AtomCloudBuilder.cpp
)

if(SEB_ENABLE_GINAC_BACKEND)
    target_sources(seb-symbolic PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic-ginac/GiNaCSymbolic.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic-ginac/RegisterSpecialFunctions.cpp
    )

    target_include_directories(seb-symbolic PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic-ginac>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/seb-symbolic-ginac>
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

target_include_directories(seb PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party
)

# Link libraries to SEB
target_link_libraries(seb PUBLIC
    seb-symbolic
    GSL::gsl
    GSL::gslcblas
)

set_target_properties(seb PROPERTIES
    EXPORT_NAME seb
)

add_library(SEB::seb ALIAS seb)

if(BUILD_DOC_TUTORIALS)
    add_subdirectory(docs/tutorials/cpp)
endif()

if(BUILD_PYTHON)
    message(STATUS "Building Python bindings")

    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    execute_process(
        COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
        OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(PYBIND11_CMAKE_DIR AND EXISTS "${PYBIND11_CMAKE_DIR}")
        set(pybind11_DIR "${PYBIND11_CMAKE_DIR}" CACHE PATH "pybind11 CMake package directory")
        list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")
    endif()
    find_package(pybind11 CONFIG REQUIRED)

    set(BINDING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/bindings/bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/bindings/bindingsTypes.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/bindings/bindingsSymbolic.cpp
    )

    pybind11_add_module(_pyseb MODULE ${BINDING_SOURCES})
    target_link_libraries(_pyseb PRIVATE SEB::seb GSL::gsl GSL::gslcblas)

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

    # scikit-build already installs into the Python package directory selected
    # by setup.py's cmake_install_dir. A second "pyseb" component here would
    # produce pyseb/pyseb/_pyseb in wheels.
    if(SKBUILD)
        set(PYSEB_PYTHON_PACKAGE_DESTINATION ".")
        set(PYSEB_PYTHON_MODEL_DESTINATION "models")
    else()
        set(PYSEB_PYTHON_PACKAGE_DESTINATION "pyseb")
        set(PYSEB_PYTHON_MODEL_DESTINATION "pyseb/models")
    endif()

    install(TARGETS _pyseb DESTINATION ${PYSEB_PYTHON_PACKAGE_DESTINATION})

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/examples/python)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pyseb/examples/python/
            DESTINATION examples/python
            FILES_MATCHING PATTERN "*.py")

    install(FILES
        ${PYSEB_BUNDLED_MODEL_FILES}
        DESTINATION ${PYSEB_PYTHON_MODEL_DESTINATION})
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

if(SEB_ENABLE_GINAC_BACKEND)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/seb-symbolic-ginac/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/seb-symbolic-ginac
            FILES_MATCHING
                PATTERN "*.hpp"
                PATTERN "build" EXCLUDE
                PATTERN "__pycache__" EXCLUDE)
endif()

install(FILES
    ${PYSEB_BUNDLED_MODEL_FILES}
    DESTINATION ${CMAKE_INSTALL_DATADIR}/pyseb/models)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/third_party/fkYAML/LICENSE.txt
        DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/pyseb
        RENAME fkYAML-LICENSE.txt)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/THIRD_PARTY_NOTICES.md
        DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/pyseb)

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
