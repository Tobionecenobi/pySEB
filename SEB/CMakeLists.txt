cmake_minimum_required(VERSION 3.18)
project(seb_tests)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required libraries
find_library(GINAC_LIBRARY NAMES ginac)
find_library(CLN_LIBRARY NAMES cln)
find_package(GSL REQUIRED)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add the test executable
add_executable(test_expression test_expression.cpp SymbolicInterface.cpp GiNaCSymbolic.cpp SpecialFunctions.cpp RegisterSpecialFunctions.cpp)

# Link libraries
target_link_libraries(test_expression
    ${GINAC_LIBRARY}
    ${CLN_LIBRARY}
    GSL::gsl
    GSL::gslcblas
)
