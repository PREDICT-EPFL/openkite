#ifndef MOBILE_ROBOT_H
#define MOBILE_ROBOT_H

#include "casadi/casadi.hpp"
#include <chrono>
#include "sys/stat.h"
#include <fstream>

enum TIRE_MODEL
{
    LINEAR,
    DUGOFF,
    BRUSH
};

struct MobileRobotProperties
{
    MobileRobotProperties()
    {
        mass = 1196.0; //1124.0;
        b    = 1.168;  //1.35;
        a    = 1.73;   //1.55;
        Iz   = mass * a * b;
        r_f  = 0.31;
        r_r  = 0.35;
        fwIz = 1.42;
        rwIz = 2.06;
        Cy_r = 240000.0;
        Cx_r = 284000.0;
        Cy_f = 226000.0;
        Cx_f = 376000.0;
        mu   = 0.8;
        rollResist = 600.0;
        rollResistSpeed = 10;
        Cxx = 1.0;
        tire_model = BRUSH;
    }
    ~MobileRobotProperties(){}

    double mass;     // mass of the car
    double b;        // distance from CoG to the rear axis
    double a;        // disance from CoG to the front axis
    double Iz;       // inertia of the car
    double r_f;      // radius of the front wheel
    double r_r;      // radius of the rear wheel
    double fwIz;     // inertia of the front wheel
    double rwIz;     // inertia of the rear wheel
    double Cy_r;     // lateral stiffness of the rear wheel
    double Cx_r;     // longitudinal stiffness of the rear wheel
    double Cy_f;     // lateral stiffness of the front wheel
    double Cx_f;     // longitudinal stiffness of the front wheel
    double mu;       // some coefficient
    double rollResist;  //rolling resistance
    double rollResistSpeed; //roll resist speed coeff
    double Cxx;      // air drag force coefficient
    TIRE_MODEL tire_model; // tire model id


};

/** Relatively simple mobile model : tricycle */
class MobileRobot
{
public:
    MobileRobot(const MobileRobotProperties &props);
    MobileRobot();
    ~MobileRobot(){}

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

#endif // MOBILE_ROBOT_H
