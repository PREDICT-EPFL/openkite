#ifndef KITE_H
#define KITE_H

#include "casadi/casadi.hpp"
#include "yaml-cpp/yaml.h"
#include "kitemath.h"
#include <chrono>
#include "sys/stat.h"
#include <fstream>

struct WindProperties {
    double WindFrom{M_PI};
    double WindSpeed{0.0};
};
struct PlaneGeometry {
    //double imuPitchOffset_deg;

    double wingSpan;
    double mac;
    double aspectRatio;
    double wingSurfaceArea;
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


    /* SidesliGENangle (beta) -> Side force (Y) */
    double CYb;

    /* SidesliGENangle (beta) -> Rolling moment (l) */
    double Cl0;
    double Clb;

    /* SidesliGENangle (beta) -> Yawing moment (n) */
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
    TetherProperties Tether;
};

struct AlgorithmProperties {
    IntType Integrator;
    double sampling_time;
};

namespace kite_utils {
/** load properties from a YAML file */
KiteProperties LoadProperties(const std::string &filename);

KiteProperties LoadMinimalProperties(const std::string &filename);


/** architecture-dependend time stamping */
typedef std::chrono::time_point<std::chrono::system_clock> time_point;

time_point get_time();

/** write a DM vector to a file */
void write_to_file(const std::string &filename, const casadi::DM &data);

/** read a vector from a file */
casadi::DM read_from_file(const std::string &filename);

/** check if file exists */
bool file_exists(const std::string &filename);

enum IdentMode {
    LONGITUDINAL,
    LATERAL,
    YAW,
    COMPLETE
};
}

class KiteDynamics {
public:
    //constructor
    KiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps,
                 const bool teth_ON = false, const bool controlsIncludeWind = false);

    KiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps,
                 const kite_utils::IdentMode &identMode,
                 const bool controlsIncludeWind = false);

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
    casadi::Function getNumericAirspeedMeas() { return this->NumAirspeedMeasured; }
    casadi::Function getNumericAeroValues() { return this->NumAeroValues; }
    casadi::Function getNumericSpecNongravForce() { return this->NumSpecNongravForce; }
    casadi::Function getNumericSpecTethForce() { return this->NumSpecTethForce; }
    casadi::Function getNumericIntegrator() { return this->NumIntegrator; }
    casadi::Function getNumericJacobian() { return this->NumJacobian; }
    casadi::Function getAeroDynamicForces() { return this->AeroDynamics; }

    /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder*/
    template<typename W, typename GEN, typename DLO, typename DLA, typename AIL, typename ELV, typename RUD>
    void getModel(GEN &g, GEN &rho,
                  W &windFrom, W &windSpeed,
                  GEN &b, GEN &c, GEN &AR, GEN &S,
                  GEN &Mass, GEN &Ixx, GEN &Iyy, GEN &Izz, GEN &Ixz,

                  GEN &e_oswald,
                  DLO &CD0,


                  DLO &CL0,
                  DLO &CLa,

                  DLO &Cm0,
                  DLO &Cma,


                  DLA &CYb,

                  GEN &Cl0,
                  DLA &Clb,

                  DLA &Cn0,
                  DLA &Cnb,


                  DLO &CLq,
                  DLO &Cmq,

                  DLA &CYp,
                  DLA &Clp,
                  DLA &Cnp,

                  DLA &CYr,
                  DLA &Clr,
                  DLA &Cnr,


                  ELV &CLde,
                  ELV &Cmde,

                  AIL &Clda,
                  AIL &Cnda,

                  RUD &CYdr,
                  RUD &Cldr,
                  RUD &Cndr,

                  casadi::SX &v, casadi::SX &w, casadi::SX &r, casadi::SX &q,
                  casadi::SX &T, casadi::SX &dE, casadi::SX &dR, casadi::SX &dA,
                  casadi::SX &v_dot, casadi::SX &w_dot, casadi::SX &r_dot, casadi::SX &q_dot,
                  casadi::SX &Va_pitot, casadi::SX &Va, casadi::SX &alpha, casadi::SX &beta,
                  casadi::SX &b_F_aero, casadi::SX &b_F_thrust, bool teth_ON, casadi::SX &b_F_tether);

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
    //numerical airspeed evaluation
    casadi::Function NumAirspeedMeasured;
    //numerical aero values
    casadi::Function NumAeroValues;
    //numerical specific nongravitational force evaluation
    casadi::Function NumSpecNongravForce;
    //numerical tether force evaluation
    casadi::Function NumSpecTethForce;
    //numerical integral evaluation
    casadi::Function NumIntegrator;
    //numerical jacobian evaluation
    casadi::Function NumJacobian;

    casadi::Function AeroDynamics;
};

class MinimalKiteDynamics {
public:

    //constructor
    MinimalKiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps,
                        const bool teth_ON, const bool controlsIncludeWind = false);

    MinimalKiteDynamics(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps,
                        const kite_utils::IdentMode &identMode,
                        const bool controlsIncludeWind = false);

    MinimalKiteDynamics() = default;

    virtual ~MinimalKiteDynamics() {}

    /** public methods */
    casadi::SX getSymbolicState() { return this->State; }
    casadi::SX getSymbolicControl() { return this->Control; }
    casadi::SX getSymbolicParameters() { return this->Parameters; }
    casadi::SX getSymbolicDynamics() { return this->SymDynamics; }
    casadi::SX getSymbolicIntegrator() { return this->SymIntegartor; }
    casadi::SX getSymbolicJacobian() { return this->SymJacobian; }

    casadi::Function getNumericDynamics() { return this->NumDynamics; }
    casadi::Function getNumericAirspeedMeas() { return this->NumAirspeedMeasured; }
    casadi::Function getNumericAeroValues() { return this->NumAeroValues; }
    casadi::Function getNumericSpecNongravForce() { return this->NumSpecNongravForce; }
    casadi::Function getNumericSpecTethForce() { return this->NumSpecTethForce; }
    casadi::Function getNumericIntegrator() { return this->NumIntegrator; }
    casadi::Function getNumericJacobian() { return this->NumJacobian; }
    casadi::Function getAeroDynamicForces() { return this->AeroDynamics; }

    /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder*/
    template<typename W, typename GEN, typename DLO, typename DLA, typename AIL, typename ELV, typename RUD>
    void getMinimalModel(GEN &g, GEN &rho,
                         W &windFrom, W &windSpeed,
                         GEN &b, GEN &c, GEN &AR, GEN &S,
                         GEN &Mass, GEN &Ixx, GEN &Iyy, GEN &Izz, GEN &Ixz,

                         DLO &e_oswald,
                         DLO &CD0,


                         DLO &CL0,
                         DLO &CLa,

                         DLO &Cm0,
                         DLO &Cma,


                         DLA &CYb,

            //GEN &Cl0,
                         DLA &Clb,

            //GEN &Cn0,
                         DLA &Cnb,


            //DLO &CLq,
            //DLO &Cmq,

            //GEN &CYp,
                         DLA &Clp,
            //DLA &Cnp,

            //DLA &CYr,
            //DLA &Clr,
                         DLA &Cnr,


            //ELV &CLde,
                         ELV &Cmde,

                         AIL &Clda,
            //GEN &Cnda,

            //RUD &CYdr,
            //RUD &Cldr,
                         RUD &Cndr,

                         casadi::SX &v, casadi::SX &w, casadi::SX &r, casadi::SX &q,
                         casadi::SX &T, casadi::SX &dE, casadi::SX &dR, casadi::SX &dA,
                         casadi::SX &v_dot, casadi::SX &w_dot, casadi::SX &r_dot, casadi::SX &q_dot,
                         casadi::SX &Va_pitot, casadi::SX &Va, casadi::SX &alpha, casadi::SX &beta,
                         casadi::SX &b_F_aero, casadi::SX &b_F_thrust, bool teth_ON, casadi::SX &b_F_tether);

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
    //numerical airspeed evaluation
    casadi::Function NumAirspeedMeasured;
    //numerical aero values
    casadi::Function NumAeroValues;
    //numerical specific nongravitational force evaluation
    casadi::Function NumSpecNongravForce;
    //numerical tether force evaluation
    casadi::Function NumSpecTethForce;
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
