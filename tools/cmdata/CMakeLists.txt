cmake_minimum_required(VERSION 3.16)

project(cmdata VERSION 0.1
               DESCRIPTION "A programm to calculate contact data from gromacs trajectories for multi-eGO"
               LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMDATA cmdata)

set(CMAKE_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(${CMDATA} main.cpp)

# set build type to release if not set (check if its debug, Debug, or DEBUG)
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CMAKE_BUILD_TYPE Debug)
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
elseif(CMAKE_BUILD_TYPE MATCHES "Release" OR NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

# find gromacs
find_package(GROMACS REQUIRED NAMES gromacs gromacs_mpi gromacs_d gromacs_mpi_d HINTS "$ENV{GROMACS_DIR}")

# include source and header files
target_include_directories(${CMDATA} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

include(FetchContent)
SET(FETCHCONTENT_QUIET OFF)
SET(FETCHCONTENT_BASE_DIR ${CMAKE_CURRENT_BINARY_DIR}/cmdata_fetch)

message(STATUS "Fetching popt")
FetchContent_Declare(
  popt
  URL https://github.com/rpm-software-management/popt/archive/refs/heads/master.zip
  GIT_TAG 2bca0aa
  PATCH_COMMAND patch --directory=${FETCHCONTENT_BASE_DIR}/popt-src < ${CMAKE_CURRENT_SOURCE_DIR}/popt.patch
)
FetchContent_MakeAvailable(popt)

message(STATUS "Fetching xdrfile")
FetchContent_Declare(
  xdrfile
  URL https://github.com/multi-ego/xdrfile/archive/refs/heads/chemfiles.zip 
)
FetchContent_MakeAvailable(xdrfile)

# link libraries
target_link_libraries(${CMDATA} PRIVATE Gromacs::libgromacs xdrfile popt)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS ${CMDATA} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# build test
if(CMDATA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test/)
endif()
