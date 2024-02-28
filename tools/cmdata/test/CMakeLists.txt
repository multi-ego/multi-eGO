cmake_minimum_required(VERSION 3.16)

enable_testing()
include(GoogleTest)
# include(GoogleMock)

# compile google test

file(GLOB TEST_LIST "*/*.cpp")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

message(STATUS "Creating build for Tests ::")
foreach(TEST ${TEST_LIST})
    message(STATUS " :: creating build for ${TEST}")
    get_filename_component(T_NAME ${TEST} NAME_WE)
    add_executable(${T_NAME} ${TEST})
    target_include_directories(${T_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/src/cmdata)
    target_link_libraries(${T_NAME} GTest::gtest_main GTest::gmock_main Gromacs::libgromacs xdrfile popt)
    target_compile_features(${T_NAME} PUBLIC cxx_std_17)
    set_property(TARGET ${T_NAME} PROPERTY CXX_STANDARD 17)
    gtest_discover_tests(${T_NAME})
endforeach(TEST ${TEST_LIST})
