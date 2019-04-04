# AI accelerated PBRT
# Minimum CMake Version
cmake_minimum_required (VERSION 3.3)

# Project Name
project (AIAPBRT)

# setze language standard
set(CMAKE_CXX_STANDARD 17) 	# enable C++17 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall")

# Optimierung
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")

# Versionsnummer
set (${PROJECT_NAME}_VERSION_MAJOR 1)
set (${PROJECT_NAME}_VERSION_MINOR 0)

# hinzufuegen aller files im src Folder
file(GLOB_RECURSE src_files
	"src/*.h"
	"src/*.cpp"
)

# erstellen der Filter fuer die src_files
source_group(FILES ${src_files} TREE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# erstellen des Targets ${PROJECT_NAME}
add_executable(${PROJECT_NAME} ${src_files})
