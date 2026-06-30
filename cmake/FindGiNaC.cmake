# FindGiNaC.cmake
# Find the GiNaC library
#
# This module defines:
#  GINAC_FOUND - True if GiNaC was found
#  GINAC_INCLUDE_DIRS - The GiNaC include directories
#  GINAC_LIBRARIES - The GiNaC libraries

find_path(GINAC_INCLUDE_DIR
  NAMES ginac/ginac.h
  PATHS
    /usr/include
    /usr/local/include
)

find_library(GINAC_LIBRARY
  NAMES ginac
  PATHS
    /usr/lib
    /usr/local/lib
    /usr/lib/x86_64-linux-gnu
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GiNaC DEFAULT_MSG
  GINAC_LIBRARY GINAC_INCLUDE_DIR)

if(GINAC_FOUND)
  set(GINAC_LIBRARIES ${GINAC_LIBRARY})
  set(GINAC_INCLUDE_DIRS ${GINAC_INCLUDE_DIR})
endif()

mark_as_advanced(GINAC_INCLUDE_DIR GINAC_LIBRARY)
