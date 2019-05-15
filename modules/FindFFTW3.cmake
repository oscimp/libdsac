# Source: https://github.com/jedbrown/cmake-modules/blob/master/FindFFTW.cmake
#
# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
	# Already in cache, be silent
	set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)

# Get the core and threads library
find_library (FFTW_CORE NAMES fftw3)
find_library (FFTW_THREADS NAMES fftw3_threads)

list(INSERT FFTW_LIBRARIES 0 ${FFTW_CORE} ${FFTW_THREADS})

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
