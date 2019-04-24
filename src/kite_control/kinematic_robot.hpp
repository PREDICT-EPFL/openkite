#ifndef KINEMATIC_ROBOT_HPP
#define KINEMATIC_ROBOT_HPP

#include "casadi/casadi.hpp"
#include <chrono>
#include "sys/stat.h"
#include <fstream>
#include "mobile_robot.h"


struct KinematicRobotProperties
{
    KinematicRobotProperties(const double &_wb = 1.0) : wheel_base(_wb){}
    ~KinematicRobotProperties(){}
    double wheel_base;
};

/** Relatively simple mobile model : tricycle */
class KinematicRobot
{
public:
    KinematicRobot(const MobileRobotProperties &props);
    KinematicRobot();
    ~KinematicRobot(){}

    casadi::Function getDynamics(){return NumDynamics;}
    casadi::Function getOutputMapping(){return OutputMap;}

    casadi::Function TraceFunction;
private:
    casadi::SX state;
    casadi::SX control;
    casadi::SX Dynamics;

    casadi::Function NumDynamics;
    casadi::Function OutputMap;
};

#endif // KINEMATIC_ROBOT_HPP
