cmake_minimum_required(VERSION 2.8.3)
project(kite_model)


include_directories(include ${CASADI_INCLUDE_DIR} src/kite_math)

add_library(kitemodel kite.cpp kite.h)
target_link_libraries(kitemodel kitemath ${YAML_CPP_LIBS})

add_library(odesolver integrator.cpp integrator.h)
target_link_libraries(odesolver kitemath kitemodel)

add_executable(kite_model_test kite_model_test.cpp)
target_link_libraries(kite_model_test odesolver kitemodel ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_executable(simulator simulator.cpp simulator.h)
target_link_libraries(simulator mobile_robot kinematic_robot odesolver ${catkin_LIBRARIES})

add_dependencies(simulator openkite_generate_messages_cpp)

add_executable(homotopy_test homotopy_test.cpp)
target_link_libraries(homotopy_test kiteNMPF kitemodel odesolver)

add_executable(psarc_test psarc_test.cpp)
target_link_libraries(psarc_test kitemodel odesolver)

add_executable(psarc_simple_test psarc_simple_test.cpp)
target_link_libraries(psarc_simple_test kitemodel odesolver)
