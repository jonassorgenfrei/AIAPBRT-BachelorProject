# AI accelerated PBRT
# Minimum CMake Version
cmake_minimum_required (VERSION 3.3)

# For sanitizers
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Project Name
project (AIAPBRT)

OPTION(PBRT_FLOAT_AS_DOUBLE "Use 64-bit floats" OFF)

IF (PBRT_FLOAT_AS_DOUBLE)
  ADD_DEFINITIONS ( -D PBRT_FLOAT_AS_DOUBLE )
ENDIF()

OPTION(PBRT_SAMPLED_SPECTRUM "Use SampledSpectrum rather than RGBSpectrum" OFF)

IF (PBRT_SAMPLED_SPECTRUM)
  ADD_DEFINITIONS ( -D PBRT_SAMPLED_SPECTRUM )
ENDIF()

if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/ext/glog/cmake")
  message(FATAL_ERROR "The glog submodule directory is missing! "
    "You probably did not clone the project with --recursive, or you first checked out "
    "pbrt before it was added. It is possible to recover by running "
    "\"git submodule update --init --recursive\"")
endif()

# setze language standard
set(CMAKE_CXX_STANDARD 17) 	# enable C++17 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall")

# Optimierung
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")

# Versionsnummer
set (${PROJECT_NAME}_VERSION_MAJOR 1)
set (${PROJECT_NAME}_VERSION_MINOR 0)

###########################################################################
# glog

SET(WITH_GFLAGS OFF CACHE BOOL "Use gflags")
SET(BUILD_SHARED_LIBS OFF CACHE BOOL " " FORCE)
IF(WIN32)
  ADD_DEFINITIONS( -D GOOGLE_GLOG_DLL_DECL= )
ENDIF()
ADD_SUBDIRECTORY(src/ext/glog)
SET_PROPERTY(TARGET glog logging_unittest demangle_unittest utilities_unittest stl_logging_unittest PROPERTY FOLDER "ext")
INCLUDE_DIRECTORIES (
  src/ext/glog/src
  ${CMAKE_BINARY_DIR}/src/ext/glog
)

###########################################################################

SET(ALL_PBRT_LIBS
 glog
)

# hinzufuegen aller files im src Folder
file(GLOB_RECURSE src_files
	"src/*.h"
	"src/*.cpp"
)

# erstellen der Filter fuer die src_files
source_group(FILES ${src_files} TREE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# erstellen des Targets ${PROJECT_NAME}
add_executable(${PROJECT_NAME} ${src_files})

TARGET_LINK_LIBRARIES (${PROJECT_NAME} ${ALL_PBRT_LIBS} )