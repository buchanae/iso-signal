# - Try to find Imagebamtools
# Once done, this will define
#
#  bamtools_FOUND - system has bamtools
#  bamtools_INCLUDE_DIRS - the bamtools include directories
#  bamtools_LIBRARIES - link these to use bamtools

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(bamtools_PKGCONF bamtools)

# Include dir
find_path(bamtools_INCLUDE_DIR
  NAMES api/BamAux.h
  PATHS ${bamtools_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES bamtools
)

# Finally the library itself
find_library(bamtools_LIBRARY
  NAMES bamtools libbamtools
  PATHS ${bamtools_PKGCONF_LIBRARY_DIRS}
  PATH_SUFFIXES bamtools
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(bamtools_PROCESS_INCLUDES bamtools_INCLUDE_DIR bamtools_INCLUDE_DIRS)
set(bamtools_PROCESS_LIBS bamtools_LIBRARY bamtools_LIBRARIES)
libfind_process(bamtools)
