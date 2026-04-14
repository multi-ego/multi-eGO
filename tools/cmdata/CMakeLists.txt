cmake_minimum_required(VERSION 3.16)

project(cmdata VERSION 0.1
               DESCRIPTION "A programm to calculate contact data from gromacs trajectories for multi-eGO"
               LANGUAGES CXX C)

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

# check if gromacs version is at least 2026
string(REPLACE "." ";" GMX_VER_LIST "${GROMACS_VERSION}")
list(GET GMX_VER_LIST 0 GMX_MAJOR_VER)

if(GMX_MAJOR_VER STRGREATER_EQUAL "2026")
    message("Using GROMACS version greater-equal than 2026")
    add_compile_definitions(GMXVGE2026)
else()
    message("Using GROMACS version less than 2026")
endif()

find_package(HDF5 COMPONENTS C CXX)
if(HDF5_FOUND)
  message(STATUS "HDF5 found, enabling HDF5 support")
  add_definitions(-DUSE_HDF5)
  target_include_directories(${CMDATA} PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${CMDATA} PRIVATE ${HDF5_CXX_LIBRARIES})
endif()

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
  PATCH_COMMAND patch --directory=${FETCHCONTENT_BASE_DIR}/popt-src -p0 < ${CMAKE_CURRENT_SOURCE_DIR}/popt.patch
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(popt)

message(STATUS "Fetching xdrfile")
FetchContent_Declare(
  xdrfile
  URL https://github.com/multi-ego/xdrfile/archive/refs/heads/chemfiles.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# xdrfile's CMakeLists.txt calls project() without VERSION; set the policy
# default so the warning is suppressed in that subdirectory.
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
FetchContent_MakeAvailable(xdrfile)

# link libraries
target_link_libraries(${CMDATA} PRIVATE Gromacs::libgromacs xdrfile popt)

# Workaround for macOS 26+ SDK: libm.tbd was removed because math functions are
# now fully part of libSystem. GROMACS (and other deps) may export an absolute
# path to libm.tbd in their CMake interface, which causes "No rule to make target"
# errors at build time. Strip those references from every imported target after
# all find_package / FetchContent calls have run.
if(APPLE)
  macro(strip_libm _tgt)
    if(TARGET ${_tgt})
      foreach(_prop INTERFACE_LINK_LIBRARIES LINK_LIBRARIES)
        get_target_property(_libs ${_tgt} ${_prop})
        if(_libs)
          list(FILTER _libs EXCLUDE REGEX ".*libm(\\.tbd)?$")
          set_property(TARGET ${_tgt} PROPERTY ${_prop} "${_libs}")
        endif()
      endforeach()
    endif()
  endmacro()
  strip_libm(Gromacs::libgromacs)
  strip_libm(xdrfile)
  strip_libm(popt)
  strip_libm(${CMDATA})
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS ${CMDATA} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# build test
option(CMDATA_BUILD_TESTS "Build and register unit tests (no GROMACS required)" OFF)
if(CMDATA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test/)
endif()
