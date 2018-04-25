#ifndef KITEEKF_H
#define KITEEKF_H

#include "kite.h"
#include "chrono"

#define DEPRECATED

class KiteEKF
{
public:
    KiteEKF(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps); DEPRECATED
    KiteEKF(std::shared_ptr<KiteDynamics> obj_Kite);                                DEPRECATED
    KiteEKF(const casadi::Function &_Dynamics, const casadi::Function &_Jacobian);

    virtual ~KiteEKF(){}

    static const casadi::DM SIGMA_Q;
    static const casadi::DM SIGMA_R;
    static const casadi::DM SIGMA_V;
    static const casadi::DM SIGMA_W;

    static const casadi::DM DEFAULT_PROCESS_COVARIANCE;
    static const casadi::DM DEFAULT_MEASUREMENT_COVARIANCE;
    static const casadi::DM DEFAULT_MEASUREMENT_MATRIX;

    void setProcessCovariance(const casadi::DM &_W){this->W = _W;}
    void setMeasurementCovariance(const casadi::DM &_V){this->V = _V;}
    void setEstimationCovariance(const casadi::DM &_P){this->Cov_Est = _P;}
    void setEstimation(const casadi::DM &_estimation){this->State_Est = _estimation;}
    void setControl(const casadi::DM &_control){this->Control = _control;}
    void setTime(const double &_current_time){this->tstamp = _current_time;}

    casadi::DM getEstimation(){return this->State_Est;}
    casadi::DM getEstimationCovariance(){return this->Cov_Est;}
    double getTimeStamp(){return this->tstamp;}

    void estimate(const casadi::DM &measurement, const double &_tstamp);
    void _estimate(const casadi::DM &measurement, const double &_dt);

    /** @brief Kalman filter equations */
    void propagate(const double & _dt);

private:
    std::shared_ptr<KiteDynamics> Kite;
    casadi::DM State_Est;
    casadi::DM Cov_Est;
    casadi::DM Control;

    casadi::Function m_Integrator;
    casadi::Function m_Jacobian;

    casadi::DM W;
    casadi::DM V;
    casadi::DM H;
    double tstamp;
};


#endif // KITEEKF_H
