cmake_minimum_required(VERSION 3.27)
project(SCPD_Project)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(include/algorithms)
include_directories(include/data_structures)
include_directories(include/utilities)
include_directories(tests)

add_executable(SCPD_Project
        include/algorithms/Algorithm1.h
        include/algorithms/Algorithm2.h
        include/data_structures/nbody.h
        include/utilities/Timer.h
        include/utilities/UtilityFunctions.h
        include/utilities/vec3.h
        src/algorithms/Parallel/Algorithm1_Parallel.cc
        src/algorithms/Parallel/Algorithm2_Parallel.cc
        src/algorithms/Sequential/Algorithm2_Sequential.cc
        src/algorithms/Sequential/all_pairs_sequential.cc
        src/data_structures/DataStructure.cc
        src/utilities/Timer.cc
        src/utilities/UtilityFunctions.cc
        src/main.cc
        tests/algorithms/Algorithm1_Test.cpp
        tests/algorithms/Algorithm2_Test.cpp
        tests/test_utils.h
        tests/data_structures/nbody_test.cpp
        include/data_structures/particle.h
        src/data_structures/nbody.cc
)
