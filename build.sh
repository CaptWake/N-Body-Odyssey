#!/bin/bash
b="cmake-build/$1"
cmake -S . -B "$b" -DCMAKE_BUILD_TYPE="$1" -DCMAKE_MAKE_PROGRAM=$(which ninja) -G Ninja
cmake --build "$b" --target SCPD_Project -j 10
