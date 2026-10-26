cmake_minimum_required(VERSION 3.25)

project(NelderMead LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(nmtest
    test.cpp
    nelder-mead.h
)

add_executable(nmtestx2y2
    test-x2-y2.cpp
    nelder-mead.h
)

