#ifndef KITENMPF_H
#define KITENMPF_H

#include "kite.h"
#include <memory>


enum CollocationType{GLOBAL, MULTIPLE_SHOOTING};

class KiteNMPF
{
public:
    //KiteNMPF(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps);
    KiteNMPF(std::shared_ptr<KiteDynamics> _Kite, const casadi::Function &_Path);

    //KiteNMPF(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps, const casadi::Function &_Path);
    virtual ~KiteNMPF(){}

    /** contsraints setters */
    void setLBX(const casadi::DM &_lbx){this->LBX = _lbx;}
    void setUBX(const casadi::DM &_ubx){this->UBX = _ubx;}

    void setLBG(const casadi::DM &_lbg){this->LBG = _lbg;}
    void setUBG(const casadi::DM &_ubg){this->UBG = _ubg;}

    void setLBU(const casadi::DM &_lbu){this->LBU = _lbu;}
    void setUBU(const casadi::DM &_ubu){this->UBU = _ubu;}

    void setStateScaling(const casadi::DM &Scaling){Scale_X = Scaling;
                                                      invSX = casadi::DM::solve(Scale_X, casadi::DM::eye(Scale_X.size1()));}
    void setControlScaling(const casadi::DM &Scaling){Scale_U = Scaling;
                                                      invSU = casadi::DM::solve(Scale_U, casadi::DM::eye(Scale_U.size1()));}

    void setReferenceVelocity(const casadi::DM &vel_ref){reference_velocity = Scale_X(14,14) * vel_ref;}

    void setPath(const casadi::SX &_path);
    void createNLP();

    void enableWarmStart(){WARM_START = true;}
    void disableWarmStart(){WARM_START = false;}
    void computeControl(const casadi::DM &_X0);
    casadi::DM findClosestPointOnPath(const casadi::DM &position, const casadi::DM &init_guess = casadi::DM(0));

    casadi::DM getOptimalControl(){return OptimalControl;}
    casadi::DM getOptimalTrajetory(){return OptimalTrajectory;}
    casadi::Function getPathFunction(){return PathFunc;}
    casadi::Function getAugDynamics(){return AugDynamics;}
    casadi::Dict getStats(){return stats;}
    bool initialized(){return _initialized;}

    double getPathError();
    double getVirtState();
    double getVelocityError();

    /** lower-upper bounds for state variables */
    static const casadi::DM DEFAULT_LBX;
    static const casadi::DM DEFAULT_UBX;

    /** lower-upper bounds for nonlinear constraints */
    static const casadi::DM DEFAULT_LBG;
    static const casadi::DM DEFAULT_UBG;

    /** lower-upper bounds for control variables */
    static const casadi::DM DEFAULT_LBU;
    static const casadi::DM DEFAULT_UBU;

private:
    std::shared_ptr<KiteDynamics> Kite;
    casadi::SX Path;
    casadi::Function PathFunc;

    casadi::SX Contraints;
    casadi::Function ContraintsFunc;

    casadi::SX reference_velocity;

    /** state box constraints */
    casadi::DM LBX, UBX;

    /** nonlinear inequality constraints */
    casadi::DM LBG, UBG;

    /** control box constraints */
    casadi::DM LBU, UBU;

    /** state and control scaling matrixces */
    casadi::DM Scale_X, invSX;
    casadi::DM Scale_U, invSU;

    /** cost function weight matrices */
    casadi::SX Q, R, W, Wq;

    casadi::DM NLP_X, NLP_LAM_G, NLP_LAM_X;
    casadi::Function NLP_Solver;
    casadi::SXDict NLP;
    casadi::Dict OPTS;
    casadi::DMDict ARG;
    casadi::Dict stats;

    casadi::DM OptimalControl;
    casadi::DM OptimalTrajectory;

    unsigned NUM_SHOOTING_INTERVALS;
    bool WARM_START;
    bool _initialized;
    bool scale;

    /** TRACE FUNCTIONS */
    casadi::Function DynamicsFunc;
    casadi::Function DynamicConstraints;
    casadi::Function PerformanceIndex;
    casadi::Function CostFunction;
    casadi::Function PathError;
    casadi::Function VelError;

    casadi::Function AugJacobian;
    casadi::Function AugDynamics;
};

#endif // KITENMPF_H
