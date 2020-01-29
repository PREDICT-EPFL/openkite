#ifndef KITE_H
#define KITE_H

#include "casadi/casadi.hpp"
#include "yaml-cpp/yaml.h"
#include "kitemath.h"
#include <chrono>
#include "sys/stat.h"
#include <fstream>

struct WindProperties {
    double WindFrom_deg{180.0};
    double WindSpeed{0.0};
};

struct PlaneGeometry {
    double imuPitchOffset_deg;

    double wingSpan;
    double mac;
    double aspectRatio;
    double wingSurfaceArea;
//    double TaperRatio;
//    double HTailsurface;
//    double TailLeverArm;
//    double FinSurfaceArea;
//    double FinLeverArm;
//    double AerodynamicCenter;
};

struct PlaneInertia {
    double mass;
    double Ixx;
    double Iyy;
    double Izz;
    double Ixz;
};

struct PlaneAerodynamics {
    double e_oswald;

    /* Static Drag force (D) (total Drag force of aircraft) */
    double CD0;


    /* Angle of attack (alpha) -> Lift force (L) (total Lift force of aircraft) */
    double CL0;
    double CLa;

    /* Angle of attack (alpha) -> Pitching moment (m) */
    double Cm0;
    double Cma;


    /* Sideslip angle (beta) -> Side force (Y) */
    double CYb;

    /* Sideslip angle (beta) -> Rolling moment (l) */
    double Cl0;
    double Clb;

    /* Sideslip angle (beta) -> Yawing moment (n) */
    double Cn0;
    double Cnb;


    /* Pitchrate (q) -> Lift force (L), Pitching moment (m) */
    double CLq;
    double Cmq;

    /* Rollrate (p) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
    double CYp;
    double Clp;
    double Cnp;

    /* Yawrate (r) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
    double CYr;
    double Clr;
    double Cnr;


    /** Aerodynamic effects of control **/
    /* Elevator deflection (de) -> Lift force (L), Pitching moment (m), Drag force (D) */
    double CLde;
    double Cmde;

    /* Aileron deflection (da) -> Rolling moment (l), Yawing moment (n) */
    double Clda;
    double Cnda;

    /* Rudder deflection (dr) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
    double CYdr;
    double Cldr;
    double Cndr;
};

//struct PlaneAerodynamics {
//    double e_oswald;
//
//    /* Static Drag force (D) (total Drag force of aircraft) */
//    double CD0;
//
//    struct AoA {
//        /* Angle of attack (alpha) -> Lift force (L) (total Lift force of aircraft) */
//        double CL0;
//        double CLa;
//
//        /* Angle of attack (alpha) -> Pitching moment (m) */
//        double Cm0;
//        double Cma;
//    } aoa;
//
//
//    struct Sideslip {
//        /* Sideslip angle (beta) -> Side force (Y) */
//        double CYb;
//
//        /* Sideslip angle (beta) -> Rolling moment (l) */
//        double Cl0;
//        double Clb;
//
//        /* Sideslip angle (beta) -> Yawing moment (n) */
//        double Cn0;
//        double Cnb;
//    } ss;
//
//
//    struct Rates {
//
//        struct Pitchrate {
//            /* Pitchrate (q) -> Lift force (L), Pitching moment (m) */
//            double CLq;
//            double Cmq;
//        } pitch;
//
//        struct Rollrate {
//            /* Rollrate (p) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
//            double CYp;
//            double Clp;
//            double Cnp;
//        } roll;
//
//        struct Yawrate {
//            /* Yawrate (r) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
//            double CYr;
//            double Clr;
//            double Cnr;
//        } yaw;
//
//    } rates;
//
//
//    struct Controls {
//        /** Aerodynamic effects of control **/
//
//        struct Elevator {
//            /* Elevator deflection (de) -> Lift force (L), Pitching moment (m), Drag force (D) */
//            double CLde;
//            double Cmde;
//        } elev;
//
//        struct Ailerons {
//            /* Aileron deflection (da) -> Rolling moment (l), Yawing moment (n) */
//            double Clda;
//            double Cnda;
//        } ail;
//
//        struct Rudder {
//            /* Rudder deflection (dr) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
//            double CYdr;
//            double Cldr;
//            double Cndr;
//        } rud;
//
//    } ctrl;
//
//};
//struct Controls {
//
//    struct Elevator {
//        /* Elevator deflection (de) -> Lift force (L), Pitching moment (m), Drag force (D) */
//        double CLde;
//        double Cmde;
//    } elev;
//
//    struct Ailerons {
//        /* Aileron deflection (da) -> Rolling moment (l), Yawing moment (n) */
//        double Clda;
//        double Cnda;
//    } ail;
//
//    struct Rudder {
//        /* Rudder deflection (dr) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
//        double CYdr;
//        double Cldr;
//        double Cndr;
//    } rud;
//
//} ctrl;
//
//struct PlaneAerodynamics {
//    double e_oswald;
//
//    /* Static Drag force (D) (total Drag force of aircraft) */
//    double CD0;
//};
//
//struct PlaneAeroAoa {
//    /* Angle of attack (alpha) -> Lift force (L) (total Lift force of aircraft) */
//    double CL0;
//    double CLa;
//
//    /* Angle of attack (alpha) -> Pitching moment (m) */
//    double Cm0;
//    double Cma;
//};
//
//struct PlaneAeroSs {
//    /* Sideslip angle (beta) -> Side force (Y) */
//    double CYb;
//
//    /* Sideslip angle (beta) -> Rolling moment (l) */
//    double Cl0;
//    double Clb;
//
//    /* Sideslip angle (beta) -> Yawing moment (n) */
//    double Cn0;
//    double Cnb;
//};
//
//
//struct PlaneAeroRatePitch {
//    /* Pitchrate (q) -> Lift force (L), Pitching moment (m) */
//    double CLq;
//    double Cmq;
//};
//
//struct PlaneAeroRateRoll {
//    /* Rollrate (p) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
//    double CYp;
//    double Clp;
//    double Cnp;
//};
//
//struct PlaneAeroRateYaw {
//    /* Yawrate (r) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
//    double CYr;
//    double Clr;
//    double Cnr;
//};
//
///** Aerodynamic effects of control **/
//struct PlaneAeroControlElev {
//    /* Elevator deflection (de) -> Lift force (L), Pitching moment (m), Drag force (D) */
//    double CLde;
//    double Cmde;
//};
//
//struct PlaneAeroControlAilerons {
//    /* Aileron deflection (da) -> Rolling moment (l), Yawing moment (n) */
//    double Clda;
//    double Cnda;
//};
//
//struct PlaneAeroControlRudder {
//    /* Rudder deflection (dr) -> Side force (Y), Rolling moment (l), Yawing moment (n) */
//    double CYdr;
//    double Cldr;
//    double Cndr;
//};

struct TetherProperties {
    double length;
    double Ks;
    double Kd;
    double rx;
    double ry;
    double rz;
};

struct KiteProperties {
    std::string name;
    WindProperties Wind;
    PlaneGeometry Geometry;
    PlaneInertia Inertia;

    PlaneAerodynamics Aerodynamics;

//    PlaneAerodynamics aero;
//    PlaneAeroAoa aero_aoa;
//    PlaneAeroSs aero_ss;
//    PlaneAeroRatePitch aero_rate_pitch;
//    PlaneAeroRateRoll aero_rate_roll;
//    PlaneAeroRateYaw aero_rate_yaw;
//    PlaneAeroControlElev aero_ctrl_elev;
//    PlaneAeroControlAilerons aero_ctrl_ail;
//    PlaneAeroControlRudder aero_ctrl_rud;

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

    enum IdentMode {
        LONGITUDINAL = 1,
        LATERAL = 2
    };

    //constructor
    KiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps);

    KiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps, const IdentMode &identMode);

    KiteDynamics() = default;

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

    template<typename P, typename LO, typename LA>
    void getModel(P &g, P &rho,
                  P &windFrom_deg, P &windSpeed,
                  P &b, P &c, P &AR, P &S,
                  P &Mass, P &Ixx, P &Iyy, P &Izz, P &Ixz,

                  LO &imuPitchOffset_deg,


                  P &e_o,
                  LO &CD0,


                  LO &CL0,
                  LO &CLa,

                  LO &Cm0,
                  LO &Cma,


                  LA &CYb,

                  P &Cl0,
                  LA &Clb,

                  P &Cn0,
                  LA &Cnb,


                  LO &CLq,
                  LO &Cmq,

                  LA &CYp,
                  LA &Clp,
                  LA &Cnp,

                  LA &CYr,
                  LA &Clr,
                  LA &Cnr,


                  LO &CLde,
                  LO &Cmde,

                  LA &Clda,
                  LA &Cnda,

                  P &CYdr,
                  P &Cldr,
                  P &Cndr,

                  casadi::SX &v, casadi::SX &w, casadi::SX &r, casadi::SX &q,
                  casadi::SX &T, casadi::SX &dE, casadi::SX &dR, casadi::SX &dA,
                  casadi::SX &v_dot, casadi::SX &w_dot, casadi::SX &r_dot, casadi::SX &q_dot,
                  casadi::SX &Faero_b, casadi::SX &T_b);

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
