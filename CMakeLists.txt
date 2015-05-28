cmake_minimum_required(VERSION 2.8.4)
project(RSA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(SOURCE_FILES main.cpp BigNum.cpp BigNum.h RSA.cpp RSA.h)
add_executable(RSA ${SOURCE_FILES} BigNum.cpp BigNum.h RSA.cpp RSA.h)