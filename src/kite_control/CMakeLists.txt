cmake_minimum_required(VERSION 2.8.3)
project(kite_control)

include_directories(${CMAKE_SOURCE_DIR}/src/kite_model ${CMAKE_SOURCE_DIR}/src/kite_estimation)

add_library(kiteNMPF kiteNMPF.cpp kiteNMPF.h)
target_link_libraries(kiteNMPF kitemodel)

# mobile robot model
add_library(mobile_robot mobile_robot.cpp mobile_robot.h)
target_link_libraries(mobile_robot kitemodel)

add_library(kinematic_robot kinematic_robot.cpp kinematic_robot.hpp)
target_link_libraries(kinematic_robot kitemodel)

#modile robot test
add_executable(car_model_test car_model_test.cpp)
target_link_libraries(car_model_test mobile_robot kinematic_robot kitemodel odesolver)

# stochastic kite model
#add_library(stochastic_kite stochastic_kite.cpp stochastic_kite.hpp)
#target_link_libraries(stochastic_kite kitemodel)

add_executable(kite_control_test kite_control_test.cpp)
target_link_libraries(kite_control_test kiteNMPF kiteEKF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_executable(kite_identification_test kite_identification_test.cpp)
target_link_libraries(kite_identification_test kiteNMPF)

#add_executable(nmpf_node nmpf_node.cpp nmpf_node.hpp)
#target_link_libraries(nmpf_node kiteNMPF odesolver mobile_robot ${catkin_LIBRARIES})

add_executable(nmpc_node nmpc_node.cpp nmpc_node.hpp)
target_link_libraries(nmpc_node odesolver mobile_robot kinematic_robot ${catkin_LIBRARIES})

#add_executable(spline_following spline_following.cpp )
#target_link_libraries(spline_following kiteNMPF)

add_executable(manifold_learning manifold_learning.cpp)
target_link_libraries(manifold_learning ${CASADI_LIBRARIES})

#add_dependencies(nmpf_node openkite_generate_messages_cpp)
