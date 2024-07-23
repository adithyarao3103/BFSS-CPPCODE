cmake_minimum_required (VERSION 2.6)

project (BFSS_SIMULATIONS)
set (PROJECT_VERSION_MAJOR 0)

set (PROJECT_VERSION_MINOR 1)


list(APPEND CMAKE_CONFIGURATION_TYPES Profile)
list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
     "Add the configurations that we need"
     FORCE)
set(CMAKE_CXX_FLAGS_PROFILE "-g -pg" CACHE STRING
    "Flags used by the C++ compiler during maintainer builds."
    FORCE )
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG_MODE -Wall -Weffc++ -Wno-unknown-pragmas")

option(USE_OPENMP "If this is activated the flags for openmp are automatically detected and
       openmp is used." ON)

set(HAVE_OPENMP OFF)

if(USE_OPENMP)
find_package(OpenMP)
if(OPENMP_FOUND)
 set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
 set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
 set(HAVE_OPENMP ON)
endif(OPENMP_FOUND)
endif(USE_OPENMP)

configure_file(
${BFSS_SIMULATIONS_SOURCE_DIR}/config.h.cmake
${BFSS_SIMULATIONS_BINARY_DIR}/config.h)

INCLUDE_DIRECTORIES(${BFSS_SIMULATIONS_BINARY_DIR})

ADD_SUBDIRECTORY (src)
