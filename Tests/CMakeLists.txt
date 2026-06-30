# Minimum CMake version
cmake_minimum_required(VERSION 3.18)

# Project name
project(pyseb_tests)

# Enable testing
enable_testing()

# Find required packages
find_package(GSL REQUIRED)
find_package(GTest REQUIRED)

# Use pkg-config to find GiNaC
include(FindPkgConfig)
pkg_check_modules(GINAC REQUIRED ginac)
pkg_check_modules(CLN REQUIRED cln)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/../SEB)
include_directories(${CMAKE_SOURCE_DIR}/../src)
include_directories(${GINAC_INCLUDE_DIRS})
include_directories(${CLN_INCLUDE_DIRS})

# Ensure GiNaC implementation code paths are enabled in headers
add_compile_definitions(USE_GINAC USE_GINAC_IMPL)

# Locate libseb and libseb-symbolic from the same build directory.
set(SEB_BUILD_DIR "" CACHE PATH "Build directory containing libseb and libseb-symbolic")
set(SEB_BUILD_HINTS
    ${CMAKE_SOURCE_DIR}/../build-ninja-ginac
    ${CMAKE_SOURCE_DIR}/../build
    ${CMAKE_SOURCE_DIR}/../build-cmake
)

if(SEB_BUILD_DIR)
    list(PREPEND SEB_BUILD_HINTS ${SEB_BUILD_DIR})
endif()

foreach(candidate_dir IN LISTS SEB_BUILD_HINTS)
    unset(CANDIDATE_SEB_LIBRARY CACHE)
    unset(CANDIDATE_SEB_SYMBOLIC_LIBRARY CACHE)

    find_library(CANDIDATE_SEB_LIBRARY
        NAMES seb libseb.a
        HINTS ${candidate_dir}
        NO_DEFAULT_PATH
    )

    find_library(CANDIDATE_SEB_SYMBOLIC_LIBRARY
        NAMES seb-symbolic libseb-symbolic.a
        HINTS ${candidate_dir}
        NO_DEFAULT_PATH
    )

    if(CANDIDATE_SEB_LIBRARY AND CANDIDATE_SEB_SYMBOLIC_LIBRARY)
        set(SEB_LIBRARY ${CANDIDATE_SEB_LIBRARY})
        set(SEB_SYMBOLIC_LIBRARY ${CANDIDATE_SEB_SYMBOLIC_LIBRARY})
        break()
    endif()
endforeach()

if(NOT SEB_LIBRARY OR NOT SEB_SYMBOLIC_LIBRARY)
    message(FATAL_ERROR
        "Could not find libseb and libseb-symbolic. Build the core library first (e.g. in ../build-ninja-ginac or ../build)."
    )
endif()

get_filename_component(SEB_LIBRARY_DIR ${SEB_LIBRARY} DIRECTORY)

# Link directories
link_directories(${SEB_LIBRARY_DIR})
link_directories(${GINAC_LIBRARY_DIRS})
link_directories(${CLN_LIBRARY_DIRS})

# Define test executables
add_executable(utility_test utility_test.cpp)
add_executable(structure_test structure_tests.cpp)
add_executable(validation_test validation_tests.cpp)
add_executable(seb_symbolic_test seb_symbolic_test.cpp)

# Add sample test executables
add_executable(sample_solid_cylinder_test Sample_SolidCylinder_test.cpp)
add_executable(sample_solid_sphere_test Sample_SolidSphere_test.cpp)
add_executable(sample_solid_spherical_shell_test Sample_SolidSphericalShell_test.cpp)
add_executable(sample_thin_disk_test Sample_ThinDisk_test.cpp)

# Add SymbolInterface implementation
list(APPEND SOURCE_FILES ${CMAKE_SOURCE_DIR}/../SEB/SymbolInterface.cpp)

# Link the libraries
target_link_libraries(utility_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

target_link_libraries(structure_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

target_link_libraries(validation_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

target_link_libraries(seb_symbolic_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

target_link_libraries(sample_solid_cylinder_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

target_link_libraries(sample_solid_sphere_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

target_link_libraries(sample_solid_spherical_shell_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

target_link_libraries(sample_thin_disk_test PRIVATE
    ${SEB_LIBRARY}
    ${SEB_SYMBOLIC_LIBRARY}
    GSL::gsl
    GSL::gslcblas
    GTest::GTest
    GTest::Main
    ${GINAC_LIBRARIES}
    ${CLN_LIBRARIES}
)

# Add compile definitions for testing
target_compile_definitions(sample_solid_cylinder_test PRIVATE TESTING)
target_compile_definitions(sample_solid_sphere_test PRIVATE TESTING)
target_compile_definitions(sample_solid_spherical_shell_test PRIVATE TESTING)
target_compile_definitions(sample_thin_disk_test PRIVATE TESTING)
target_compile_definitions(validation_test PRIVATE VALIDATION_DIR="${CMAKE_SOURCE_DIR}/../Examples/Validation/")

# Add tests
add_test(
    NAME UtilityTest
    COMMAND utility_test
)

add_test(
    NAME StructureTest
    COMMAND structure_test
)

add_test(
    NAME ValidationTest
    COMMAND validation_test
)

add_test(
    NAME SebSymbolicTest
    COMMAND seb_symbolic_test
)

# Add sample tests
add_test(
    NAME SampleSolidCylinderTest
    COMMAND sample_solid_cylinder_test
)

add_test(
    NAME SampleSolidSphereTest
    COMMAND sample_solid_sphere_test
)

add_test(
    NAME SampleSolidSphericalShellTest
    COMMAND sample_solid_spherical_shell_test
)

add_test(
    NAME SampleThinDiskTest
    COMMAND sample_thin_disk_test
)

# Set test properties
set_tests_properties(UtilityTest PROPERTIES TIMEOUT 30)
set_tests_properties(StructureTest PROPERTIES TIMEOUT 30)
set_tests_properties(ValidationTest PROPERTIES TIMEOUT 30)
set_tests_properties(SebSymbolicTest PROPERTIES TIMEOUT 30)
set_tests_properties(SampleSolidCylinderTest PROPERTIES TIMEOUT 90)
set_tests_properties(SampleSolidSphereTest PROPERTIES TIMEOUT 30)
set_tests_properties(SampleSolidSphericalShellTest PROPERTIES TIMEOUT 90)
set_tests_properties(SampleThinDiskTest PROPERTIES TIMEOUT 30)
