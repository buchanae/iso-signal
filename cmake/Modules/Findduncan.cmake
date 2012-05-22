# - Try to find Imageduncan
# Once done, this will define
#
#  duncan_FOUND - system has duncan
#  duncan_INCLUDE_DIRS - the duncan include directories
#  duncan_LIBRARIES - link these to use duncan

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(duncan_PKGCONF duncan)

# Include dir
find_path(duncan_INCLUDE_DIR
  NAMES duncan/Feature.h
  PATHS ${duncan_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(duncan_LIBRARY
  NAMES duncan libduncan
  PATHS ${duncan_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(duncan_PROCESS_INCLUDES duncan_INCLUDE_DIR duncan_INCLUDE_DIRS)
set(duncan_PROCESS_LIBS duncan_LIBRARY duncan_LIBRARIES)
libfind_process(duncan)
