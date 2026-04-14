cmake_minimum_required(VERSION 3.16)

# Test targets — intentionally have NO dependency on GROMACS or any external
# library.  Only indexing.hpp, density.hpp, and mindist.hpp are needed; all
# three are pure C++17 header-only code with no GROMACS includes.

add_executable(test_cmdata test_cmdata.cpp)
target_include_directories(test_cmdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_compile_features(test_cmdata PRIVATE cxx_std_17)
target_compile_options(test_cmdata PRIVATE -Wall -Wextra -Wno-unused-parameter)
add_test(NAME cmdata_tests COMMAND test_cmdata)
