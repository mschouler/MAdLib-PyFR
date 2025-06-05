# -------------------------------------------------------------------
# MAdLib - Copyright (C) Cenaero
#
# See the Copyright.txt and License.txt files for license information.
# -------------------------------------------------------------------

#########################################################################
# Main CMakeLists.txt for MAdLib
#########################################################################
cmake_minimum_required(VERSION 3.18)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()


# -- PLATFORM --
# --------------

#enable_language(CXX)

# Additional folders to search for CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Platforms/")

if(DEFINED ENV{PLATFORM})
	if (EXISTS "${CMAKE_SOURCE_DIR}/cmake/Platforms/$ENV{PLATFORM}.cmake")
	  message("Including platform $ENV{PLATFORM}")
	  include($ENV{PLATFORM})
	else()
	  message("Platform configuration file $ENV{PLATFORM}.cmake does not exist, continueing with default configuration.")
	endif()
else()
  message("PLATFORM environment variable was not defined, continueing with default configuration.")
endif()


# Flags for MSVC
# set(CMAKE_USER_MAKE_RULES_OVERRIDE
#   ${CMAKE_SOURCE_DIR}/cmake/Modules/c_flag_overrides.cmake)
# set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
#   ${CMAKE_SOURCE_DIR}/cmake/Modules/cxx_flag_overrides.cmake)

#set(CMAKE_BUILD_TYPE Release)

if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING
  "Choose build type, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel"
  )
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING
  "Choose build type, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel")
endif()


# -- DEFAULTS --
# --------------

if (NOT DEFINED MADLIB_USE_PARALLEL_DEFAULT)
   set (MADLIB_USE_PARALLEL_DEFAULT FALSE)
endif()


# -- PARALLEL SETTINGS --
# -----------------------

set(MADLIB_USE_PARALLEL ${MADLIB_USE_PARALLEL_DEFAULT} CACHE BOOL "Use sequential or parallel version - compiler and mpi library are set according to environment")

########################################################################

project(MAdLib)

# Define C++ standard
set(CMAKE_CXX_STANDARD 11)

if (MADLIB_USE_ANN)
  project(ANN)
endif()

if (MADLIB_USE_CGNS)
  project(CGNS)
endif()

if (MADLIB_USE_PARALLEL AND NOT AUTOPACK_ROOT)
  project(autopack)
endif()

##########################################################################
# Get svn revision number - stored in cmake/Modules/MAdConfig.h
##########################################################################
find_program(SVNVERSION svnversion)
mark_as_advanced(SVNVERSION)
execute_process(COMMAND ${SVNVERSION} -c ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE MADLIB_SVN_REVISION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX REPLACE ".*\:" "" MADLIB_SVN_REVISION ${MADLIB_SVN_REVISION})

#########################################################################
# MAdLib revision number
#########################################################################
set(MADLIB_VERSION_MAJOR 2)
set(MADLIB_VERSION_MINOR 4)

set(MADLIB_USE_DISTRIBUTION_DEFAULT OFF)
set(MADLIB_USE_DISTRIBUTION ${MADLIB_USE_DISTRIBUTION_DEFAULT} CACHE BOOL
"Release distribution - default: ${MADLIB_USE_DISTRIBUTION_DEFAULT}")
mark_as_advanced(MADLIB_USE_DISTRIBUTION)

if(MADLIB_USE_DISTRIBUTION)
  set(MADLIB_VERSION_PATCH 0)
else()
  set(MADLIB_VERSION_PATCH ${MADLIB_SVN_REVISION})
endif()
set(MADLIB_VERSION "${MADLIB_VERSION_MAJOR}.${MADLIB_VERSION_MINOR}.${MADLIB_VERSION_PATCH}")

#########################################################################
# Disallow in-source build
if("${MAdLib_SOURCE_DIR}"  STREQUAL "${MAdLib_BINARY_DIR}")
  message(FATAL_ERROR
    "MAdLib requires an out of source Build. Please create a separate binary directory and run CMake there.")
endif()

# -- Registers size --
set(MADLIB_USE_32BITS_ON_64BITS_PLATFORM_DEFAULT FALSE)
if(NOT WIN32)
  set(MADLIB_USE_32BITS_ON_64BITS_PLATFORM
        ${MADLIB_USE_32BITS_ON_64BITS_PLATFORM_DEFAULT}
        CACHE BOOL "Force 32 bits compilation on 64 bits platform - default:${MADLIB_USE_32BITS_ON_64BITS_PLATFORM_DEFAULT}")
else()
  set(MADLIB_USE_32BITS_ON_64BITS_PLATFORM ${MADLIB_USE_32BITS_ON_64BITS_PLATFORM_DEFAULT})
endif()

# Put REGISTERS_SIZE in cache because strange behavior if re-run cmake on auto
# compil on ernest (CMAKE_SIZEOF_VOID_P empty)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  if(MADLIB_USE_32BITS_ON_64BITS_PLATFORM)
    set(REGISTERS_SIZE 32 CACHE STRING "Registers size")
    set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS OFF)
  else()
    set(REGISTERS_SIZE 64 CACHE STRING "Registers size")
  endif()
else()
  set(REGISTERS_SIZE 64 CACHE STRING "Registers size")
endif()
mark_as_advanced(REGISTERS_SIZE)

if(MADLIB_USE_32BITS_ON_64BITS_PLATFORM)
  if(DEFINED REGISTERS_SIZE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m${REGISTERS_SIZE}")
    set(LINK_FLAGS "${LINK_FLAGS} -m${REGISTERS_SIZE}")
  endif()
endif()

# RPATH handling - execution possible from
  # use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  # the RPATH to be used when installing
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")###
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# -- Default compiler option --
# -----------------------------
if( UNIX )  # To avoid empty CMAKE_CXX_COMPILER_ID
  if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wshadow -Wno-write-strings")
  endif()
endif()

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DPROFILE")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPROFILE -DDEBUG -D DEBUG")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL XL AND
   ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "12")  # Compiler installed on AIX
  set(CMAKE_CXX_FLAGS "-q64  -qrtti=all")
  set(CMAKE_C_FLAGS "-q64  -qrtti=all")
  set(CMAKE_Fortran_FLAGS "-q64  -qrtti=all")  # 64 bits - allow dynamic_casts
  set(CMAKE_EXE_LINKER_FLAGS "-q64 -bbigtoc")
  # To take -X64 into account, must invert flags
  # See http://www.cmake.org/pipermail/cmake/2010-July/038312.html
  set(XL_STATIC_LIBRARY_FLAGS "-X64")
  set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> <LINK_FLAGS> cr <TARGET> <OBJECTS>")
  set(CMAKE_C_ARCHIVE_APPEND "<CMAKE_AR> <LINK_FLAGS> r <TARGET> <OBJECTS>")
  set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> <LINK_FLAGS> cr <TARGET> <OBJECTS>")
  set(CMAKE_CXX_ARCHIVE_APPEND "<CMAKE_AR> <LINK_FLAGS> r <TARGET> <OBJECTS>")
endif()

if(WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE) # e.g. warning for unsafe usage of fopen - proposed solution is non-standard
endif()

add_compile_options($<$<CXX_COMPILER_ID:GNU>:-Werror=vla>)

###########################################################
# Build options
###########################################################
set(MADLIB_BUILD_BENCHMARKS_DEFAULT OFF)
set(MADLIB_BUILD_BENCHMARKS ${MADLIB_BUILD_BENCHMARKS_DEFAULT} CACHE BOOL "Build all available benchmarks (default: ${MADLIB_BUILD_BENCHMARKS_DEFAULT})" )

set(MADLIB_BUILD_TESTCASES_DEFAULT ON)
set(MADLIB_BUILD_TESTCASES ${MADLIB_BUILD_TESTCASES_DEFAULT} CACHE BOOL "Build all available testcases (default: ${MADLIB_BUILD_TESTCASES_DEFAULT})" )

set(MADLIB_USE_ANN_DEFAULT OFF)
set(MADLIB_USE_ANN ${MADLIB_USE_ANN_DEFAULT} CACHE BOOL "Build/use ANN library (default: ${MADLIB_USE_ANN_DEFAULT})" )

set(MADLIB_USE_CGNS_DEFAULT OFF)
set(MADLIB_USE_CGNS ${MADLIB_USE_CGNS_DEFAULT} CACHE BOOL "Interface CGNS (default: ${MADLIB_USE_CGNS_DEFAULT})" )

set(MADLIB_USE_CGNSCPEX0045_DEFAULT OFF)
set(MADLIB_USE_CGNSCPEX0045 ${MADLIB_USE_CGNSCPEX0045_DEFAULT} CACHE BOOL "Interface CGNS with CPEX0045 (default: ${MADLIB_USE_CGNS_DEFAULT})" )

set(MADLIB_USE_GMM_DEFAULT OFF)
set(MADLIB_USE_GMM ${MADLIB_USE_GMM_DEFAULT} CACHE BOOL "Compile gmm++ linear solvers (default: ${MADLIB_USE_GMM_DEFAULT})" )

set(MADLIB_USE_MATHEX_DEFAULT ON)
set(MADLIB_USE_MATHEX ${MADLIB_USE_MATHEX_DEFAULT} CACHE BOOL "Use Mathex library (default: ${MADLIB_USE_MATHEX_DEFAULT})" )

set(MADLIB_USE_METIS_DEFAULT OFF)
set(MADLIB_USE_METIS ${MADLIB_USE_METIS_DEFAULT} CACHE BOOL "Use Metis partitionner (default: ${MADLIB_USE_METIS_DEFAULT})" )

set(MADLIB_USE_OCC_DEFAULT OFF)
set(MADLIB_USE_OCC ${MADLIB_USE_OCC_DEFAULT} CACHE BOOL "Use OpenCascade geometric model through Gmsh (default: ${MADLIB_USE_OCC_DEFAULT})" )

set(MADLIB_USE_PARMETIS_DEFAULT OFF)
set(MADLIB_USE_PARMETIS ${MADLIB_USE_PARMETIS_DEFAULT} CACHE BOOL "Use ParMetis partitionner (default: ${MADLIB_USE_PARMETIS_DEFAULT})" )

set(MADLIB_USE_PTSCOTCH_DEFAULT OFF)
set(MADLIB_USE_PTSCOTCH ${MADLIB_USE_PTSCOTCH_DEFAULT} CACHE BOOL "Use PT-SCOTCH partitionner (default: ${MADLIB_USE_PTSCOTCH_DEFAULT})" )

set(MADLIB_USE_PARSER_DEFAULT OFF)
set(MADLIB_USE_PARSER ${MADLIB_USE_PARSER} CACHE BOOL "Use Parser (restricted, CENAERO, Belgium - default: ${MADLIB_USE_PARSER})" )

set(MADLIB_USE_PETSC_DEFAULT OFF)
set(MADLIB_USE_PETSC ${MADLIB_USE_PETSC_DEFAULT} CACHE BOOL "Use PETSc (default: ${MADLIB_USE_PETSC_DEFAULT})" )

set(MADLIB_USE_SPARSKIT_DEFAULT OFF)
set(MADLIB_USE_SPARSKIT ${MADLIB_USE_SPARSKIT_DEFAULT} CACHE BOOL "Use Sparskit linear solver (default: ${MADLIB_USE_SPARSKIT_DEFAULT})" )

if (NOT DEFINED MADLIB_BUILD_SHARED_LIBRARY_DEFAULT)
  set(MADLIB_BUILD_SHARED_LIBRARY_DEFAULT ON)
endif()
set(MADLIB_BUILD_SHARED_LIBRARY ${MADLIB_BUILD_SHARED_LIBRARY_DEFAULT} CACHE BOOL
  "Build MAdlib shared library (default:${MADLIB_BUILD_SHARED_LIBRARY_DEFAULT})")

if(NOT DEFINED MADLIB_BUILD_STATIC_LIBRARY_DEFAULT)
  set(MADLIB_BUILD_STATIC_LIBRARY_DEFAULT OFF)
endif()
set(MADLIB_BUILD_STATIC_LIBRARY ${MADLIB_BUILD_STATIC_LIBRARY_DEFAULT} CACHE BOOL
  "Build MAdlib static library (default:${MADLIB_BUILD_STATIC_LIBRARY_DEFAULT})")

#########################################################################
# Install options
#########################################################################

set(MADLIB_INSTALL_EXTENDED_API_DEFAULT OFF)
set(MADLIB_INSTALL_EXTENDED_API ${MADLIB_INSTALL_EXTENDED_API_DEFAULT} CACHE BOOL "Install all header files in the API instead of the 'official' ones (default: ${MADLIB_INSTALL_EXTENDED_API_DEFAULT}) ")

if( EXISTS "${CMAKE_INSTALL_PREFIX}/include_all" )  # Security check
  install(CODE "message(\"Will remove install folder '${CMAKE_INSTALL_PREFIX}'\")")
  install(CODE "file(REMOVE_RECURSE ${CMAKE_INSTALL_PREFIX})" )
endif()


#########################################################################
# Parse options
#########################################################################
set(MADLIB_SRC "")

if(NOT MADLIB_BUILD_SHARED_LIBRARY AND NOT MADLIB_BUILD_STATIC_LIBRARY)
  message(FATAL_ERROR "At least shared or static library have to be built:
  1) set MADLIB_BUILD_STATIC_LIBRARY to ON for the build of static lib.
  2) set MADLIB_BUILD_SHARED_LIBRARY to ON for the build of shared lib.
  3) set MADLIB_BUILD_STATIC_LIBRARY and MADLIB_BUILD_SHARED_LIBRARY to ON for the build both libraries.")
endif()

# ----- ANN -----
if(MADLIB_USE_ANN)
  add_definitions(-D_HAVE_ANN_)
  add_subdirectory (Contrib/ANN)
  include_directories(Contrib/ANN/include)
endif()

# ------ AUTOPACK ------
if(MADLIB_USE_PARALLEL)
  if(AUTOPACK_ROOT)
    find_package(AUTOPACK)
    if(AUTOPACK_FOUND)
      set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                       ${AUTOPACK_INCLUDE_DIRS})
      set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                    ${AUTOPACK_LIBRARIES})
    else()
      message("Autopack library is not found; Continue by compiling the library included within the MAdLib repository")
      add_definitions(-D_HAVE_AUTOPACK_)
      add_subdirectory (Contrib/autopack)
      include_directories(Contrib/autopack)
    endif()
  else()
    message("No AUTOPACK_ROOT path defined; Continue by compiling the library included within the MAdLib repository -- check")
    add_definitions(-D_HAVE_AUTOPACK_)
    add_subdirectory (Contrib/autopack)
    include_directories(Contrib/autopack)
  endif()
endif()

# -----BLAS/LAPACK -----

if (NOT BLAS_FOUND AND NOT LAPACK_FOUND)
  if( MSVC OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND DEFINED ENV{MKLROOT} )
    set(BLA_VENDOR "Intel10_64lp_seq")
    set(MADLIB_USE_MKL TRUE CACHE INTERNAL "")
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "SunOS")
    set(BLA_VENDOR "SunPerf")
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static" OR ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-dynamic")
    set(BLA_VENDOR "IBMESSL")
  elseif(APPLE)
    set(BLA_VENDOR "Apple")
  endif()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

if(BLAS_FOUND AND LAPACK_FOUND)
  set(MADLIB_USE_BLAS_LAPACK TRUE)
  set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static" 
      OR ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-dynamic")
    find_package(Threads)
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  endif()
  if(MADLIB_USE_MKL)
    add_definitions(-D_HAVE_MKL_)                                 
    set(BLASLAPACK_INCLUDE_DIRS "${MKL_ROOT}/include")
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                    ${BLASLAPACK_INCLUDE_DIRS})
  endif()
endif()

# ----- CGNS -----

set(HDF5_DIR_DEFAULT $ENV{HDF5DIR})
set(HDF5_DIR ${HDF5_DIR_DEFAULT} CACHE PATH "path to hdf5 installation")
find_library(HDF5_LIB hdf5 PATHS ${HDF5_DIR} PATH_SUFFIXES lib NO_DEFAULT_PATH)
find_library(HDF5_LIB hdf5 PATHS ENV HDF5DIR PATH_SUFFIXES lib NO_DEFAULT_PATH)
find_library(HDF5_LIB hdf5)
mark_as_advanced(HDF5_DIR)
mark_as_advanced(HDF5_LIB)

set(CGNS_DIR_DEFAULT $ENV{CGNSDIR})
set(CGNS_DIR ${CGNS_DIR_DEFAULT} CACHE PATH "path to cgns installation")
find_library(CGNS_LIB cgns PATHS ${CGNS_DIR} PATH_SUFFIXES lib NO_DEFAULT_PATH)
find_library(CGNS_LIB cgns PATHS ENV CGNSDIR PATH_SUFFIXES lib NO_DEFAULT_PATH)
find_library(CGNS_LIB cgns)
mark_as_advanced(CGNS_DIR)
mark_as_advanced(CGNS_LIB)

find_path(CGNS_INC "cgnslib.h" PATHS ${CGNS_DIR} PATH_SUFFIXES include NO_DEFAULT_PATH)
find_path(CGNS_INC "cgnslib.h" PATHS ENV CGNSDIR PATH_SUFFIXES include NO_DEFAULT_PATH)
find_path(CGNS_INC "cgnslib.h")
mark_as_advanced(CGNS_INC)

if (MADLIB_USE_CGNS)
  if (CGNS_LIB AND CGNS_INC AND HDF5_LIB)
    set(CGNS_SRC "")
    set(CGNS_API "")
    set(CGNSDIR "")
    add_subdirectory(CGNS)
    set(MADLIB_HAVE_CGNS TRUE)
    list(APPEND MADLIB_EXTERNAL_LIBRARIES ${HDF5_LIB})
    list(APPEND MADLIB_EXTERNAL_LIBRARIES ${CGNS_LIB})
		list(APPEND MADLIB_API_FILES ${CGNS_API})
		include_directories(${CGNS_INC})
		include_directories(${CGNSDIR})

    add_definitions(-D_HAVE_CGNS_)
    if(MADLIB_USE_CGNSCPEX0045)
      add_definitions(-D_HAVE_CPEX0045_)
    endif()
  else()
    if (NOT CGNS_LIB)
      message("CGNS library not found")
    endif()
    if (NOT HDF5_LIB)
      message("HDF5 library not found")
    endif()
    if (NOT CGNS_INC)
      message("CGNS includes not found")
    endif()
    message("Continueing without CGNS (MADLIB_HAVE_CGNS set to 'FALSE')")
    set(MADLIB_USE_CGNS FALSE)
  endif()
endif()

# ------ MATHEX ------
if(MADLIB_USE_MATHEX)
  add_definitions(-D_HAVE_MATHEX_)
  set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                   ${CMAKE_SOURCE_DIR}/Contrib/mathex)
endif()

# ------ METIS ------
if(MADLIB_USE_METIS)
  find_package(METIS)
  if(METIS_FOUND)
    add_definitions(-D_HAVE_METIS_)
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                     ${METIS_INCLUDE_DIRS})
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                  ${METIS_LIBRARIES})
  else()
    message("Metis library is not found; Continue with MADLIB_USE_METIS set to 'FALSE'")
    set(MADLIB_USE_METIS FALSE)
  endif()
endif()

# ------ OCC -----
if(MADLIB_USE_OCC)
  find_package(OCC)
  if(OCC_FOUND)
    add_definitions(-D_HAVE_OCC_)
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                     ${OCC_INCLUDE_DIR})
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                  ${OCC_LIBRARIES})
  else()
    message("OpenCasCade library is not found; Continue with MADLIB_USE_OCC set to 'FALSE'")
    set(MADLIB_USE_OCC FALSE)
  endif()
endif()

# ------ PARALLEL ------
if(MADLIB_USE_PARALLEL)

	find_package(MPI REQUIRED)

  if(MPI_FOUND)

    add_definitions(-D_HAVE_MPI_)
    add_definitions(-DPARALLEL)
    add_definitions(-DMPICH_IGNORE_CXX_SEEK)
    add_definitions(-DMPIPP_H)

		include_directories(SYSTEM ${MPI_INCLUDE_PATH})

		link_libraries(${MPI_LIBRARIES})

		set(COMPILE_FLAGS ${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(LINK_FLAGS ${LINK_FLAGS} ${MPI_LINK_FLAGS})

    if(AUTOPACK_FOUND)
      set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                    ${AUTOPACK_LIBRARIES})
    endif()

  else()
    message("MPI library is not found")
    set(MADLIB_USE_PARALLEL FALSE)
  endif()
endif()

# ------ PARMETIS -----
if(MADLIB_USE_PARMETIS AND MADLIB_USE_PARALLEL)
  find_package(PARMETIS)
  if(PARMETIS_FOUND)
    add_definitions(-D_HAVE_PARMETIS_)
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                     ${PARMETIS_INCLUDE_DIRS})
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                  ${PARMETIS_LIBRARIES} ${METIS_LIBRARIES})
  else()
    message("ParMetis library is not found; Continue with MADLIB_USE_PARMETIS set to 'FALSE'")
    set(MADLIB_USE_PARMETIS FALSE)
  endif()
endif()

# ------ PTSCOTCH -----
if(MADLIB_USE_PTSCOTCH AND MADLIB_USE_PARALLEL)
  find_package(PTSCOTCH)
  if(PTSCOTCH_FOUND)
    add_definitions(-D_HAVE_PTSCOTCH_)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -lrt -pthread")
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                     ${PTSCOTCH_INCLUDE_DIRS})
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                  ${PTSCOTCH_LIBRARIES})
  else()
    message("PT-SCOTCH library is not found; Continue with MADLIB_USE_PTSCOTCH set to 'FALSE'")
    set(MADLIB_USE_PTSCOTCH FALSE)
  endif()
endif()

if(MADLIB_USE_PARMETIS OR MADLIB_USE_PTSCOTCH)
  add_definitions(-D_HAVE_BALANCING_)
endif()

# ----- PARSER ------
if(MADLIB_USE_PARSER)
  find_package(PARSER)
  if(PARSER_FOUND)
    add_definitions(-D_HAVE_PARSER_)
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                     ${PARSER_INCLUDE_DIRS})
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                  ${PARSER_LIBRARIES})
  else()
    message("Parser library is not found; Continue with MADLIB_USE_PARSER set to 'FALSE'")
    set(MADLIB_USE_PARSER FALSE)
  endif()
endif()

# ----- PETSC -----
if(MADLIB_USE_PETSC)
  find_package(PETSC)
  if(PETSC_FOUND)
    add_definitions(-D_HAVE_PETSC_)
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                     ${PETSC_INCLUDE_DIRS})
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                  ${PETSC_LIBRARIES})
  else()
    message("PETSC library is not found; Continue with MADLIB_USE_PETSC set to 'FALSE'")
    set(MADLIB_USE_PETSC FALSE)
  endif()
endif()

# ----- SPARSKIT -----
if(MADLIB_USE_SPARSKIT)
  find_package(SPARSKIT)
  if(SPARSKIT_FOUND)
    add_definitions(-D_HAVE_SPARSKIT_)
    set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                     ${SPARSKIT_INCLUDE_DIRS})
    set(MADLIB_EXTERNAL_LIBRARIES ${MADLIB_EXTERNAL_LIBRARIES}
                                  ${SPARSKIT_LIBRARIES})
  else()
    message("Sparskit library is not found; Continue with MADLIB_USE_SPARSKIT set to 'FALSE'")
    set(MADLIB_USE_SPARSKIT FALSE)
  endif()
endif()

# ----- BOOST Headers -----
find_package(Boost 1.61.0)
if(Boost_FOUND)
  set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                    ${Boost_INCLUDE_DIRS})
  set(HAVE_BOOST true)
elseif(NOT HAVE_BOOST)
  # If Boost is not installed, we download the version 1.70.0
  message(STATUS "Boost is not found; Continue with included Boost RTree (version 1.70.0)")
  configure_file(cmake/Modules/EP_Boost.cmake EP_Boost/CMakeLists.txt)
  # Build the boost sub-project
  execute_process(COMMAND ${CMAKE_COMMAND} -G ${CMAKE_GENERATOR} .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/EP_Boost
  )
  # Download and extract library
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/EP_Boost
  )
  set(HAVE_BOOST true)
endif()

# ----- VTK -----
# Set the path to the VTK folder
set(VTK_DIR "/path/to/vtk/build")

# Find VTK package
find_package(VTK)
if(VTK_FOUND)
  set(MADLIB_EXTERNAL_INCLUDE_DIRS ${MADLIB_EXTERNAL_INCLUDE_DIRS}
                                    ${VTK_LIBRARIES})
  message(STATUS "VTK is found")
endif()


#########################################################################
# Build sub-libraries
#########################################################################

set(SUB_DIRECTORIES ${SUB_DIRECTORIES} Mesh Model Common Adapt Transformation Mapping CInterface Solution SearchTree LevelSet Localisation)

if(MADLIB_USE_PARALLEL AND (MADLIB_USE_PTSCOTCH OR MADLIB_USE_PARMETIS))
  set(SUB_DIRECTORIES ${SUB_DIRECTORIES} Balancing)
endif()

foreach(DIR ${SUB_DIRECTORIES})
  set(${DIR}_SRC "")
  add_subdirectory(${DIR})
  list(APPEND MADLIB_SRC ${${DIR}_SRC})
  list(APPEND MADLIB_INCLUDE_DIRS ${DIR} ${${DIR}_INCLUDE_DIRS})
endforeach()

if(MSVC)
  foreach(DIR ${MADLIB_INCLUDE_DIRS})
    file(GLOB MADLIB_DIR_HEADERS ${DIR}/*.h)
    list(APPEND MADLIB_HEADERS ${MADLIB_DIR_HEADERS})
  endforeach()

  set(MADLIB_SRC ${MADLIB_SRC} ${MADLIB_HEADERS})

else()

  set(MADLIB_SRC ${MADLIB_SRC})

endif()

set(MADLIB_INCLUDE_DIRS ${MADLIB_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/cmake/Modules)


#########################################################################
# Build final library
#########################################################################
include_directories(${MAdLib_BINARY_DIR}) #to be seen later
include_directories(${MADLIB_EXTERNAL_INCLUDE_DIRS})
include_directories(${MADLIB_INCLUDE_DIRS})
if(MADLIB_BUILD_STATIC_LIBRARY)
  set(MADLIB_EXTERNAL_LIBRARIES_STATIC ${MADLIB_EXTERNAL_LIBRARIES})

  if (MADLIB_USE_ANN)
    add_library(ANNSTATIC STATIC ${ANN_SRC})
    set_target_properties(ANNSTATIC PROPERTIES OUTPUT_NAME "ANN")
    set(MADLIB_EXTERNAL_LIBRARIES_STATIC ${MADLIB_EXTERNAL_LIBRARIES_STATIC} "ANNSTATIC")
    if(WIN32)
      set_target_properties(ANNSTATIC PROPERTIES PREFIX "lib")
    endif()
  endif()

  if (MADLIB_USE_PARALLEL AND NOT AUTOPACK_FOUND)
    add_library(autopackSTATIC STATIC ${AUTOPACK_SRC})
    set_target_properties(autopackSTATIC PROPERTIES OUTPUT_NAME "autopack")
    set(MADLIB_EXTERNAL_LIBRARIES_STATIC ${MADLIB_EXTERNAL_LIBRARIES_STATIC} "autopackSTATIC")
    if(WIN32)
      set_target_properties(autopackSTATIC PROPERTIES PREFIX "lib")
    endif()
  endif()

	if (MADLIB_USE_CGNS)
	 	add_library(MAdLibCGNSSTATIC STATIC ${CGNS_SRC})
    set_target_properties(MAdLibCGNSSTATIC PROPERTIES OUTPUT_NAME "MAdLibCGNS")
		set(MADLIB_MODULE_LIBRARIES_STATIC ${MADLIB_MODULE_LIBRARIES_STATIC} "MAdLibCGNSSTATIC")
    if(WIN32)
      set_target_properties(MAdLibCGNSSTATIC PROPERTIES PREFIX "lib")
    endif()

    # add the include dirs needed for boost
    if(Boost_FOUND)
      target_include_directories(MAdLibCGNSSTATIC SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
    else()
      target_include_directories(MAdLibCGNSSTATIC SYSTEM PUBLIC Contrib/EP_Boost/src)
    endif()
  endif()
 
  add_library(MAdLibSTATIC STATIC ${MADLIB_SRC})
  set_target_properties(MAdLibSTATIC PROPERTIES OUTPUT_NAME "MAdLib")
  if(WIN32)
    set_target_properties(MAdLibSTATIC PROPERTIES PREFIX "lib")
    set(_MADLIB_STATIC_LINKING_ 1)
  endif()
   
  # add the include dirs needed for boost
  if(Boost_FOUND)
    target_include_directories(MAdLibSTATIC SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
  else()
    target_include_directories(MAdLibSTATIC SYSTEM PUBLIC Contrib/EP_Boost/src)
  endif()

endif()

if(MADLIB_BUILD_SHARED_LIBRARY)
  set(MADLIB_EXTERNAL_LIBRARIES_SHARED ${MADLIB_EXTERNAL_LIBRARIES})

  if (MADLIB_USE_ANN)
    add_library(ANN SHARED ${ANN_SRC})
    set(MADLIB_EXTERNAL_LIBRARIES_SHARED ${MADLIB_EXTERNAL_LIBRARIES_SHARED} "ANN")
  endif()

  if (MADLIB_USE_PARALLEL AND NOT AUTOPACK_FOUND)
    add_library(autopack SHARED ${AUTOPACK_SRC})
    set(MADLIB_EXTERNAL_LIBRARIES_SHARED ${MADLIB_EXTERNAL_LIBRARIES_SHARED} "autopack")
  endif()

  if (MADLIB_USE_CGNS)
     add_library(MAdLibCGNS SHARED ${CGNS_SRC})
     set(MADLIB_MODULE_LIBRARIES_SHARED ${MADLIB_MODULE_LIBRARIES_SHARED} "MAdLibCGNS")
     # add the include dirs needed for boost
     if(Boost_FOUND)
       target_include_directories(MAdLibCGNS SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
     else()
       target_include_directories(MAdLibCGNS SYSTEM PUBLIC Contrib/EP_Boost/src)
     endif()
  endif()

  add_library(MAdLib SHARED ${MADLIB_SRC})

  if(MSVC)
    # add the include dirs and external libraries, needed for VS projects
    target_include_directories(MAdLib PUBLIC "${MADLIB_EXTERNAL_INCLUDE_DIRS}")
    target_link_libraries(MAdLib "${MADLIB_EXTERNAL_LIBRARIES}")
  endif()
  
  
  # add the include dirs needed for boost
  if(Boost_FOUND)
    target_include_directories(MAdLib SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
  else()
    target_include_directories(MAdLib SYSTEM PUBLIC Contrib/EP_Boost/src)
  endif()

endif()

if(BLAS_LINKER_FLAGS} AND {LAPACK_LINKER_FLAGS})
  if(MADLIB_BUILD_STATIC_LIBRARY)
    set_target_properties(MAdLibSTATIC PROPERTIES LINK_FLAGS ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
  endif()
  if(MADLIB_BUILD_SHARED_LIBRARY)
    set_target_properties(MAdLib PROPERTIES LINK_FLAGS ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
  endif()
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL XL)  # Compiler installed on AIX
  if(MADLIB_BUILD_STATIC_LIBRARY)
    set_target_properties(MAdLibSTATIC PROPERTIES STATIC_LIBRARY_FLAGS "${XL_STATIC_LIBRARY_FLAGS}")
  endif()
endif()

if (VTK_FOUND)
    target_include_directories(MAdLib SYSTEM PUBLIC ${VTK_INCLUDE_DIRS})
    target_link_libraries(MAdLib ${VTK_LIBRARIES})
endif()

#########################################################################
# Testing
#########################################################################
add_subdirectory(cmake/DashboardScripts)

#########################################################################
# Packaging
#########################################################################
set(MADLIB_PACKAGE_SHORT_NAME "MAdLib")
set(MADLIB_PACKAGE_BUGREPORT "contrib@madlib.be")
set(PACKAGE_FULL_NAME  "${PACKAGE_SHORT_NAME}")
set(MADLIB_PACKAGE_VERSION "${MADLIB_VERSION}")
set(MADLIB_PACKAGE_STRING "${PACKAGE_FULL_NAME} ${PACKAGE_VERSION}")
set(MADLIB_PACKAGE_TARNAME  "${PACKAGE_FULL_NAME}")
set(VERSION ${MADLIB_VERSION})
set(GENERIC_DLL_VERSION ${VERSION} )
set(MADLIB_SOVERSION "${MADLIB_VERSION_MAJOR}.${MADLIB_VERSION_MINOR}")
if(MADLIB_BUILD_SHARED_LIBRARY)
  set_target_properties(MAdLib PROPERTIES SOVERSION ${MADLIB_SOVERSION} VERSION
    ${GENERIC_DLL_VERSION})
endif()
if(MADLIB_BUILD_STATIC_LIBRARY)
  set_target_properties(MAdLibSTATIC PROPERTIES SOVERSION ${MADLIB_SOVERSION}
    VERSION ${GENERIC_DLL_VERSION})
endif()

#########################################################################
# Install files
#########################################################################
# Default install path
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install CACHE PATH "Install path prefix,
prepended onto install directories." FORCE)
endif()

# --- Library ---
if(MADLIB_BUILD_SHARED_LIBRARY)
  if(WIN32)
    install(TARGETS MAdLib ARCHIVE DESTINATION bin)  # Install .dll in bin folder on Windows
    install(TARGETS MAdLib RUNTIME DESTINATION lib)  # Install .lib in lib folder on Windows
    if (MADLIB_USE_CGNS)
      install(TARGETS MAdLibCGNS ARCHIVE DESTINATION bin)
      install(TARGETS MAdLibCGNS RUNTIME DESTINATION lib)
    endif()
    if (MADLIB_USE_ANN)
      install(TARGETS ANN ARCHIVE DESTINATION bin/Contrib)
      install(TARGETS ANN RUNTIME DESTINATION lib/Contrib)
    endif()
    if (MADLIB_USE_MATHEX)
      install(TARGETS MATHEX ARCHIVE DESTINATION bin/Contrib)
      install(TARGETS MATHEX RUNTIME DESTINATION lib/Contrib)
    endif()
    if (MADLIB_USE_PARALLEL AND NOT AUTOPACK_FOUND)
      install(TARGETS autopack ARCHIVE DESTINATION bin/Contrib)
      install(TARGETS autopack RUNTIME DESTINATION lib/Contrib)
    endif()
  else()
    install(TARGETS MAdLib LIBRARY DESTINATION lib)  # Install .so in lib folder on Linux
    if (MADLIB_USE_CGNS)
      install(TARGETS MAdLibCGNS LIBRARY DESTINATION lib)
    endif()
    if (MADLIB_USE_ANN)
      install(TARGETS ANN LIBRARY DESTINATION lib/Contrib)
    endif()
    if (MADLIB_USE_PARALLEL AND NOT AUTOPACK_FOUND)
      install(TARGETS autopack LIBRARY DESTINATION lib/Contrib)
    endif()
  endif()
endif()

if(MADLIB_BUILD_STATIC_LIBRARY)
  install(TARGETS MAdLibSTATIC ARCHIVE DESTINATION lib)  # Install .a (Linux)/static .lib (Windows) in lib folder
  if (MADLIB_USE_CGNS)
    install(TARGETS MAdLibCGNSSTATIC ARCHIVE DESTINATION lib)
  endif()
  if (MADLIB_USE_ANN)
    install(TARGETS ANNSTATIC ARCHIVE DESTINATION lib/Contrib)
  endif()
  if (MADLIB_USE_PARALLEL AND NOT AUTOPACK_FOUND)
    install(TARGETS autopackSTATIC ARCHIVE DESTINATION lib/Contrib)
  endif()
endif()

# --- API ---
set(MADLIB_API_FILES ${MADLIB_API_FILES} ${MAdLib_BINARY_DIR}/cmake/Modules/MAdConfig.h)

if (MADLIB_USE_CGNS)
	 set(MADLIB_API_FILES ${MADLIB_API_FILES} ${CGNS_API})
endif()

install(FILES  ${MADLIB_API_FILES} DESTINATION include/MAdLib)

if (MADLIB_USE_ANN)
  install (FILES ${ANN_HEADERS} DESTINATION include/Contrib/ANN)
endif()

if (MADLIB_USE_MATHEX)
  install(FILES  ${CMAKE_SOURCE_DIR}/Contrib/mathex/mathex.h DESTINATION include/Contrib/mathex)
endif()

if (MADLIB_USE_PARALLEL AND NOT AUTOPACK_FOUND)
  install(FILES  ${CMAKE_SOURCE_DIR}/Contrib/autopack/autopack.h DESTINATION include/Contrib/autopack)
endif()

# --- Extended API ---

if(MADLIB_INSTALL_EXTENDED_API)
  install(FILES  ${MADLIB_EXTENDED_API_FILES} ${MADLIB_API_FILES} DESTINATION include_all/MAdLib)
endif()

# --- Documentation ---
install(FILES  ${CMAKE_SOURCE_DIR}/README
               ${CMAKE_SOURCE_DIR}/License.txt
               ${CMAKE_SOURCE_DIR}/Copyright.txt
               ${CMAKE_SOURCE_DIR}/Credits.txt
        DESTINATION ${CMAKE_INSTALL_PREFIX})

#########################################################################
# doxygen Documentation
#########################################################################
add_subdirectory(Documentation)

#########################################################################
# Benchmarks and Testcases
#########################################################################
set(CTEST_USE_LAUNCHERS 0)
include(CTest)  # To use CTest input files - must be in main CMakeLists.txt
mark_as_advanced(BUILD_TESTING)
mark_as_advanced(DART_TESTING_TIMEOUT)

configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake.in
               ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

add_subdirectory(Benchmarks)
add_subdirectory(Testcases)


#########################################################################
# Configure file
#########################################################################
CONFIGURE_FILE(${MAdLib_SOURCE_DIR}/cmake/Modules/MAdConfig.h.cmake ${MAdLib_BINARY_DIR}/cmake/Modules/MAdConfig.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/cmake/Modules)

configure_file(${CMAKE_SOURCE_DIR}/Copyright.txt.in
               ${CMAKE_SOURCE_DIR}/Copyright.txt @ONLY)

#########################################################################
# Print to screen the final options list
#########################################################################

# Options
list(APPEND OPTIONS_LIST ANN BLAS_LAPACK CGNS GMM MATHEX METIS OCC PARALLEL PARMETIS PTSCOTCH PARSER PETSC SPARSKIT STATIC)

foreach(OPTION ${OPTIONS_LIST})
  if(MADLIB_USE_${OPTION})
    list(APPEND OPTIONS_TRUE ${OPTION})
  else()
    list(APPEND OPTIONS_FALSE ${OPTION})
  endif()

  if(MADLIB_USE_${OPTION})
    if(NOT MADLIB_USE_${OPTION}_DEFAULT)
      list(APPEND NOT_DEFAULT_OPTION ${OPTION})
    endif()
  else()
    if(MADLIB_USE_${OPTION}_DEFAULT)
      list(APPEND NOT_DEFAULT_OPTION ${OPTION})
    endif()
  endif()
endforeach()

# Builds
list(APPEND BUILDS_LIST BENCHMARKS DOXYGEN TESTCASES)

foreach(BUILD ${BUILDS_LIST})
  if(MADLIB_BUILD_${BUILD})
    list(APPEND BUILDS_ON ${BUILD})
  endif()
endforeach()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS_ALL "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_EXE_LINKER_FLAGS_ALL
      "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(CMAKE_CXX_FLAGS_ALL
      "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  set(CMAKE_EXE_LINKER_FLAGS_ALL
      "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS_ALL
      "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_EXE_LINKER_FLAGS_ALL
      "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
else()
  set(CMAKE_CXX_FLAGS_ALL "Unknown")
  set(CMAKE_EXE_LINKER_FLAGS_ALL "Unknown")
endif()

#find_package(OpenMP)
#if(OPENMP_FOUND)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

get_property(MADLIB_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS)

message("")
message("--------------------------------------------------------------------")
message("=> Following options are used: ")
message("- TRUE: ${OPTIONS_TRUE}")
message("- FALSE: ${OPTIONS_FALSE}")
message("")
message("=> Non default options: ${NOT_DEFAULT_OPTION}")
message("")
message("=> Build(s) in 'all' target: MAdLib ${BUILDS_ON}")
message("")
message("=> Compilation:")
message("   - CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message("   - C++:")
message("     - Compiler: ${CMAKE_CXX_COMPILER}")
if(MADLIB_USE_PARALLEL)
  message("     - MPI Compiler: ${MPI_CXX_COMPILER}")
endif()
message("     - Flags: ${CMAKE_CXX_FLAGS_ALL}")
message("   - Link:")
message("     - Flags for exe: ${CMAKE_EXE_LINKER_FLAGS_ALL}")
message("   - Definitions: ${MADLIB_DEFINITIONS}")
message("--------------------------------------------------------------------")
message("")

