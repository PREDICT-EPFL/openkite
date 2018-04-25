cmake_minimum_required(VERSION 2.8.3)
project(kite_math)

include_directories(include ${CASADI_INCLUDE_DIR})

add_library(kitemath kitemath.cpp kitemath.h)
target_link_libraries(kitemath ${CASADI_LIBRARIES} )

add_subdirectory(pseudospectral)
