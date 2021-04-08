cmake_minimum_required(VERSION 2.8.3)
project(kite_model)


include_directories(include ${CASADI_INCLUDE_DIR} src/kite_math)

add_library(kitemodel kite.cpp kite.h)
target_link_libraries(kitemodel kitemath ${YAML_CPP_LIBRARY})

add_library(odesolver integrator.cpp integrator.h)
target_link_libraries(odesolver kitemath kitemodel)

#add_executable(kite_model_test kite_model_test.cpp)
#target_link_libraries(kite_model_test odesolver kitemodel ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

#add_executable(simulator simulator.cpp simulator.h)
#target_link_libraries(simulator kitemodel odesolver ${catkin_LIBRARIES})

#add_dependencies(simulator openkite_generate_messages_cpp)

