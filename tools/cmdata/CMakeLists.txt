cmake_minimum_required(VERSION 3.16)

project(cmdata VERSION 0.1
               DESCRIPTION "A programm to calculate contact data from gromacs trajectories for multi-eGO"
               LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMDATA cmdata)

add_executable(${CMDATA} main.cpp)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

# find gromacs
find_package(GROMACS REQUIRED NAMES gromacs gromacs_mpi gromacs_d gromacs_mpi_d HINTS "$ENV{GROMACS_DIR}" )
target_link_libraries(${CMDATA} PRIVATE Gromacs::libgromacs)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})

# find gmxapi
find_package(gmxapi REQUIRED CONFIG PATHS "$ENV{GROMACS_DIR}" )

# include source and header files
target_include_directories(${CMDATA} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# add subdirectories
#add_subdirectory(src/cmdata)

# add subdirectories
add_subdirectory(src/popt)
target_include_directories(${CMDATA} PRIVATE src/popt)

#
add_subdirectory(src/xdrfile/)
target_include_directories(${CMDATA} PRIVATE src/xdrfile/include)

target_link_libraries(${CMDATA} PRIVATE xdrfile Gromacs::gmxapi popt)