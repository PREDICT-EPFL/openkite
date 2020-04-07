cmake_minimum_required(VERSION 2.8.3)
project(openkite)

set(CMAKE_MODULE_PATH /usr/local/share/cmake/Modules/)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake)


find_package(catkin REQUIRED COMPONENTS roscpp rospy std_msgs message_generation
                                        geometry_msgs sensor_msgs roslint tf2_ros
                                        visualization_msgs serial)
find_package(VRPN REQUIRED)
find_package(CASADI REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(YAML_CPP REQUIRED yaml-cpp>=0.5.2)

if(APPLE)
  ##Due to warnings raised in CasADi
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O3")
  set(YAML_CPP_LIBS yaml-cpp)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O3")
  set(YAML_CPP_LIBS yaml-cpp)
endif()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")
set(CMAKE_VERBOSE_MAKEFILE “ON”)

set(CMAKE_BUILD_TYPE "Debug")

#detect 32bit system
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "Target architecture 64 bits")
else()
    message(STATUS "Target architecture 32 bits")
endif()


#if(ARM_ENABLE_VECTORIZATION AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm" )
#    message(STATUS "EIGEN: ARM VECTORIZATION ENABLED")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpu=neon")
#endif()

## System dependencies are found with CMake's conventions
find_package(Boost REQUIRED COMPONENTS system signals)

## Generate messages in the 'msg' folder
add_message_files(DIRECTORY msg FILES aircraft_controls.msg
                                      mpc_diagnostic.msg)

## Generate services in the 'srv' folder
add_service_files(DIRECTORY srv FILES simple.srv )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
generate_messages(DEPENDENCIES std_msgs sensor_msgs geometry_msgs)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
   CATKIN_DEPENDS message_runtime
   INCLUDE_DIRS src/kite_model src/kite_control src/kite_estimation src/nodes
#  LIBRARIES openkite
#  CATKIN_DEPENDS other_catkin_pkg
#  DEPENDS system_lib
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations

include_directories(src/kite_model src/nodes src/kite_control src/kite_estimation src/kite_math
                    ${catkin_INCLUDE_DIRS} ${VRPN_INCLUDE_DIR} ${CASADI_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR})


## Declare a C++ executable
add_executable(casadi_test src/casadi_test.cpp)
target_link_libraries(casadi_test ${CASADI_LIBRARIES} ${YAML_CPP_LIBS} ${catkin_LIBRARIES})
#target_link_libraries(lidar_driver ${catkin_LIBRARIES})

add_subdirectory(src/kite_model)
add_subdirectory(src/kite_control)
add_subdirectory(src/kite_estimation)
add_subdirectory(src/kite_math)
add_subdirectory(src/nodes)

#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
install(PROGRAMS
   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

## Mark executables and/or libraries for installation
install(TARGETS casadi_test
   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
 )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_openkite.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)

