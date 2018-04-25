cmake_minimum_required(VERSION 2.8.3)
project(kite_estimation)

include_directories(${CMAKE_SOURCE_DIR}/src/kite_model)

add_library(kiteEKF kiteEKF.cpp kiteEKF.h)
target_link_libraries(kiteEKF kitemodel)

add_executable(ekf_node ekf_node.cpp ekf_node.h)
target_link_libraries(ekf_node kiteEKF ${catkin_LIBRARIES})

add_dependencies(ekf_node openkite_generate_messages_cpp)
