cmake_minimum_required(VERSION 3.27)
project(SCPD_Project)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include_directories(include)
include_directories(include/algorithms)
include_directories(include/data_structures)
include_directories(include/utilities)
include_directories(tests)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
        argparse
        GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

add_executable(SCPD_Project
        include/data_structures/nbody.h
        include/utilities/vec3.h
        src/main.cc
        tests/test_utils.h
        include/simulation.h
        include/utilities/data_format.h
        include/data_structures/sequential_ap.h
        src/data_structures/sequential_ap.cc
        include/utilities/time_utils.h
)

target_link_libraries(SCPD_Project argparse)

add_executable(
        NBody_test
        tests/data_structures/nbody_test.cpp
        include/utilities/vec3.h
        include/data_structures/nbody.h
        include/data_structures/sequential_ap.h
        src/data_structures/sequential_ap.cc
)


target_link_libraries(
        NBody_test
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(NBody_test)
