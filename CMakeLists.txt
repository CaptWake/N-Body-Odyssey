cmake_minimum_required(VERSION 3.27)
project(SCPD_Project)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
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

add_executable(SCPD_Project
        include/algorithms/Algorithm1.h
        include/algorithms/Algorithm2.h
        include/data_structures/nbody.h
        include/utilities/vec3.h
        src/algorithms/Parallel/Algorithm1_Parallel.cc
        src/algorithms/Parallel/Algorithm2_Parallel.cc
        src/algorithms/Sequential/Algorithm2_Sequential.cc
        src/algorithms/Sequential/all_pairs_sequential.cc
        src/main.cc
        tests/algorithms/Algorithm1_Test.cpp
        tests/algorithms/Algorithm2_Test.cpp
        tests/test_utils.h
        include/data_structures/body.h
        src/data_structures/nbody.cc
        include/simulation.h
        include/utilities/data_format.h
)

add_executable(
        NBody_test
        tests/data_structures/nbody_test.cpp
        include/utilities/vec3.h
        include/data_structures/nbody.h
        include/data_structures/body.h
        src/data_structures/nbody.cc
)

target_link_libraries(
        NBody_test
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(NBody_test)
