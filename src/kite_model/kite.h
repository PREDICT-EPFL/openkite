#ifndef KITE_H
#define KITE_H

#include "casadi/casadi.hpp"
#include "yaml-cpp/yaml.h"
#include "kitemath.h"
#include <chrono>
#include "sys/stat.h"
#include <fstream>

struct WindProperties {
    double WindFrom_deg;
    double WindSpeed;
};

struct PlaneGeometry {
    double ImuPitchOffset_deg;

    double WingSpan;
    double MAC;
    double AspectRatio;
    double WingSurfaceArea;
    double TaperRatio;
    double HTailsurface;
    double TailLeverArm;
    double FinSurfaceArea;
    double FinLeverArm;
    double AerodynamicCenter;
};

struct PlaneInertia {
    double Mass;
    double Ixx;
    double Iyy;
    double Izz;
    double Ixz;
};

struct PlaneAerodynamics {
    double CL0;
    double CL0_tail;
    double CLa_total;
    double CLa_wing;
    double CLa_tail;
    double e_oswald;

    double CD0_total;
    double CD0_wing;
    double CD0_tail;
    double CYb;
    double CYb_vtail;
    double Cm0;
    double Cma;
    double Cn0;
    double Cnb;
    double Cl0;
    double Clb;

    double CLq;
    double Cmq;
    double CYr;
    double Cnr;
    double Clr;
    double CYp;
    double Clp;
    double Cnp;

    double CLde;
    double CYdr;
    double Cmde;
    double Cndr;
    double Cldr;
    double CDde;
    double Clda;
    double Cnda;
};

struct TetherProperties {
    double length;
    double Ks;
    double Kd;
    double rx;
    double ry;
    double rz;
};

struct KiteProperties {
    std::string Name;
    WindProperties Wind;
    PlaneGeometry Geometry;
    PlaneInertia Inertia;
    PlaneAerodynamics Aerodynamics;
    TetherProperties Tether;
};

struct AlgorithmProperties {
    IntType Integrator;
    double sampling_time;
};


namespace kite_utils {
    /** load properties from a YAML file */
    KiteProperties LoadProperties(const std::string &filename);

    /** architecture-dependend time stamping */
    typedef std::chrono::time_point<std::chrono::system_clock> time_point;

    time_point get_time();

    /** write a DM vector to a file */
    void write_to_file(const std::string &filename, const casadi::DM &data);

    /** read a vector from a file */
    casadi::DM read_from_file(const std::string &filename);

    /** check if file exists */
    bool file_exists(const std::string &filename);

}

class KiteDynamics {
public:
    //constructor
    KiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps);

    KiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps, const bool &id);

    virtual ~KiteDynamics() {}

    /** public methods */
    casadi::SX getSymbolicState() { return this->State; }

    casadi::SX getSymbolicControl() { return this->Control; }

    casadi::SX getSymbolicParameters() { return this->Parameters; }

    casadi::SX getSymbolicDynamics() { return this->SymDynamics; }

    casadi::SX getSymbolicIntegrator() { return this->SymIntegartor; }

    casadi::SX getSymbolicJacobian() { return this->SymJacobian; }

    casadi::Function getNumericDynamics() { return this->NumDynamics; }

    casadi::Function getNumericNumSpecNongravForce() { return this->NumSpecNongravForce; }

    casadi::Function getNumericIntegrator() { return this->NumIntegrator; }

    casadi::Function getNumericJacobian() { return this->NumJacobian; }

    casadi::Function getAeroDynamicForces() { return this->AeroDynamics; }

private:
    //state variables
    casadi::SX State;
    //control variables
    casadi::SX Control;
    //parameters
    casadi::SX Parameters;
    //symbolic equations of motion
    casadi::SX SymDynamics;
    //symbolic equations of integartor
    casadi::SX SymIntegartor;
    //symbolic expression for system jacobian
    casadi::SX SymJacobian;

    //numerical dynamics evaluation
    casadi::Function NumDynamics;
    //numerical specific nongravitational force evaluation
    casadi::Function NumSpecNongravForce;
    //numerical integral evaluation
    casadi::Function NumIntegrator;
    //numerical jacobian evaluation
    casadi::Function NumJacobian;

    casadi::Function AeroDynamics;
};

/** 6-DoF Kinematics of a Rigid Body */
class RigidBodyKinematics {
public:
    RigidBodyKinematics(const AlgorithmProperties &AlgoProps);

    virtual ~RigidBodyKinematics() {}

    casadi::Function getNumericIntegrator() { return NumIntegartor; }

    casadi::Function getNumericJacobian() { return NumJacobian; }

    casadi::Function getNumericDynamcis() { return NumDynamics; }

private:
    casadi::SX state;
    AlgorithmProperties algo_props;

    /** numerical integrator */
    casadi::Function NumIntegartor;
    /** numerical evaluation of Jacobian */
    casadi::Function NumJacobian;
    /** numerical evaluation of system dynamics */
    casadi::Function NumDynamics;
};


#endif // KITE_H
