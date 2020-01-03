cmake_minimum_required(VERSION 2.8.3)
project(nodes)

add_library(simple_vrpn_client simple_vrpn_client.cpp)
target_link_libraries(simple_vrpn_client ${catkin_LIBRARIES} ${VRPN_LIBRARIES})

add_executable(simple_logger simple_logger.cpp)
target_link_libraries(simple_logger ${catkin_LIBRARIES})

add_executable(optitrack_client optitrack_client.cpp)
target_link_libraries(optitrack_client simple_vrpn_client)

add_executable(control_proxy_node control_proxy_node.cpp control_proxy_node.hpp)
target_link_libraries(control_proxy_node ${catkin_LIBRARIES})

message(WARNING "${catkin_LIBRARIES}")

add_executable(kite_visualization_node kite_visualization_node.cpp)
target_link_libraries(kite_visualization_node kitemath ${catkin_LIBRARIES})

add_executable(transport_delay transport_delay.cpp)
target_link_libraries(transport_delay ${catkin_LIBRARIES})

add_executable(transport_delay_test transport_delay_test.cpp)
target_link_libraries(transport_delay_test ${catkin_LIBRARIES})


add_dependencies(optitrack_client openkite_generate_messages_cpp)
add_dependencies(control_proxy_node openkite_generate_messages_cpp)
add_dependencies(kite_visualization_node openkite_generate_messages_cpp)

install(TARGETS kite_visualization_node
        ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
        )
