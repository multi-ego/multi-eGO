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

# add popt
add_subdirectory(src/popt)
target_include_directories(${CMDATA} PRIVATE src/popt)

# add xdrfile
add_subdirectory(src/xdrfile/)
target_include_directories(${CMDATA} PRIVATE src/xdrfile/include)

# link libraries
target_link_libraries(${CMDATA} PRIVATE Gromacs::libgromacs xdrfile popt)

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS ${CMDATA} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# build test
if(CMDATA_BUILD_TESTS)
  add_subdirectory(test/)
endif()
