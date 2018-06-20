cmake_minimum_required(VERSION 2.8.3)
project(bifoiler)


include_directories(include ${CASADI_INCLUDE_DIR} src/bifoiler)

add_library(boatmath boat_math.cpp boat_math.h)
target_link_libraries(boatmath ${CASADI_LIBRARIES})

add_library(bifoiler boat_model.cpp boat_properties.cpp boat_model.h boat_properties.h)
target_link_libraries(bifoiler boatmath ${CASADI_LIBRARIES} ${YAML_CPP_LIBS})
