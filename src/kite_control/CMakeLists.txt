cmake_minimum_required(VERSION 2.8.3)
project(kite_control)

include_directories(${CMAKE_SOURCE_DIR}/src/kite_model ${CMAKE_SOURCE_DIR}/src/kite_estimation)

add_library(kiteNMPF kiteNMPF.cpp kiteNMPF.h)
target_link_libraries(kiteNMPF kitemodel)

add_library(mobile_robot mobile_robot.cpp mobile_robot.h)
target_link_libraries(mobile_robot kitemodel)

add_executable(kite_control_test kite_control_test.cpp)
target_link_libraries(kite_control_test kiteNMPF kiteEKF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_executable(kite_identification_test kite_identification_test.cpp)
target_link_libraries(kite_identification_test kiteNMPF)

add_executable(nmpf_node nmpf_node.cpp nmpf_node.hpp)
target_link_libraries(nmpf_node kiteNMPF odesolver mobile_robot ${catkin_LIBRARIES})

add_executable(spline_following spline_following.cpp )
target_link_libraries(spline_following kiteNMPF)

add_dependencies(nmpf_node openkite_generate_messages_cpp)
