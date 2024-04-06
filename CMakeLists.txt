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

# Enable AVX on AMD
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

add_executable(SCPD_Project
        include/data_structures/nbody.h
        include/utilities/vec3.h
        include/utilities/octree.h
        include/utilities/time_utils.h
        include/simulation.h
        include/data_structures/sequential_ap.h
        src/data_structures/sequential_ap.cc
        include/data_structures/sequential_ap_avx.h
        src/data_structures/sequential_ap_avx.cc
        include/data_structures/sequential_bh.h
        src/data_structures/sequential_bh.cpp
        src/main.cc
        src/utilities/octree.cpp
)

target_link_libraries(SCPD_Project argparse)

add_executable(
        NBody_test
        tests/data_structures/nbody_test.cpp
        include/utilities/vec3.h
        include/data_structures/nbody.h
        include/data_structures/sequential_ap.h
        src/data_structures/sequential_ap.cc
        include/nlohmann/json.hpp
)


target_link_libraries(
        NBody_test
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(NBody_test)
