# FindCLN.cmake
# Find the CLN library
#
# This module defines:
#  CLN_FOUND - True if CLN was found
#  CLN_INCLUDE_DIRS - The CLN include directories
#  CLN_LIBRARIES - The CLN libraries

find_path(CLN_INCLUDE_DIR
  NAMES cln/cln.h
  PATHS
    /usr/include
    /usr/local/include
)

find_library(CLN_LIBRARY
  NAMES cln
  PATHS
    /usr/lib
    /usr/local/lib
    /usr/lib/x86_64-linux-gnu
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLN DEFAULT_MSG
  CLN_LIBRARY CLN_INCLUDE_DIR)

if(CLN_FOUND)
  set(CLN_LIBRARIES ${CLN_LIBRARY})
  set(CLN_INCLUDE_DIRS ${CLN_INCLUDE_DIR})
endif()

mark_as_advanced(CLN_INCLUDE_DIR CLN_LIBRARY)
