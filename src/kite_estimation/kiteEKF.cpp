#include "kiteEKF.h"

using namespace casadi;

/** experimentaly defined values */
const DM KiteEKF::SIGMA_Q = DM::diag(DMVector{0.01, 0.05, 0.05, 0.05});
const DM KiteEKF::SIGMA_V = DM::diag(DMVector{0.5, 0.5, 0.5});
const DM KiteEKF::SIGMA_W = DM::diag(casadi::DMVector{0.5, 0.5, 0.5});
const DM KiteEKF::SIGMA_R = DM::diag(casadi::DMVector{0.5, 0.1, 0.1});

const DM KiteEKF::DEFAULT_PROCESS_COVARIANCE = pow(DM::diagcat({SIGMA_V, SIGMA_W, SIGMA_R, SIGMA_Q}), 2);
const DM KiteEKF::DEFAULT_MEASUREMENT_COVARIANCE = pow(DM::diag(DMVector{0.01,0.01,0.01, 0.0001,0.005,0.005,0.005}), 2);
const DM KiteEKF::DEFAULT_MEASUREMENT_MATRIX = DM::horzcat(DMVector{DM::zeros(7,6), DM::eye(7)});

KiteEKF::KiteEKF(const KiteProperties &KiteProps, const AlgorithmProperties &AlgoProps)
{
    /** instantiate kite object */
    Kite = std::make_shared<KiteDynamics>(KiteProps, AlgoProps);

    this->m_Integrator = Kite->getNumericIntegrator();
    this->m_Jacobian = Kite->getNumericJacobian();

    this->W = DEFAULT_PROCESS_COVARIANCE;
    this->V = DEFAULT_MEASUREMENT_COVARIANCE;
    this->H = DEFAULT_MEASUREMENT_MATRIX;
    this->Cov_Est = 10 * this->W;

    /** initialize time */
    std::chrono::time_point<std::chrono::system_clock> t_now = kite_utils::get_time();
    auto duration = t_now.time_since_epoch();
    auto microsec = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    this->tstamp = static_cast<double>(microsec) * 1e-6;
}

KiteEKF::KiteEKF(std::shared_ptr<KiteDynamics> obj_Kite)
{
    Kite = obj_Kite;

    this->m_Integrator = Kite->getNumericIntegrator();
    this->m_Jacobian = Kite->getNumericJacobian();

    this->W = DEFAULT_PROCESS_COVARIANCE;
    this->V = DEFAULT_MEASUREMENT_COVARIANCE;
    this->H = DEFAULT_MEASUREMENT_MATRIX;
    this->Cov_Est = 10 * this->W;

    /** initialize time */
    std::chrono::time_point<std::chrono::system_clock> t_now = kite_utils::get_time();
    auto duration = t_now.time_since_epoch();
    auto microsec = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    this->tstamp = static_cast<double>(microsec) * 1e-6;
}

KiteEKF::KiteEKF(const Function &_Dynamics, const Function &_Jacobian)
{
    /** instantiate model */
    this->m_Integrator = _Dynamics;
    this->m_Jacobian = _Jacobian;

    this->W = DEFAULT_PROCESS_COVARIANCE;
    this->V = DEFAULT_MEASUREMENT_COVARIANCE;
    this->H = DEFAULT_MEASUREMENT_MATRIX;
    this->Cov_Est = 10 * this->W;

    /** initialize time */
    std::chrono::time_point<std::chrono::system_clock> t_now = kite_utils::get_time();
    auto duration = t_now.time_since_epoch();
    auto microsec = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    this->tstamp = static_cast<double>(microsec) * 1e-6;

    std::cout << tstamp << "\n";
}


void KiteEKF::propagate(const double &_dt)
{
    DM new_state;
    if (m_Integrator.name().find("RK4") != std::string::npos)
    {
        DMVector out = m_Integrator(DMVector{State_Est, Control, _dt});
        new_state = out[0];
    }
    else if (m_Integrator.name().find("CVODES") != std::string::npos)
    {
        DMDict args = {{"x0", State_Est}, {"p", Control}};
        DMDict out = m_Integrator(args);
        new_state = out["xf"];
    }
    else
        std::cout << "WARNING: Unknown intergrator! \n";

    /** propagate covariance */
    DM A = m_Jacobian(DMVector{State_Est, Control})[0] * _dt + DM::eye(13);
    DM P_new = DM::mtimes(DM::mtimes(A, Cov_Est), A.T()) + W;

    State_Est = new_state;
    Cov_Est = P_new;
}

void KiteEKF::estimate(const DM &measurement, const double &_tstamp)
{
    /** prediction step */
    double dt = _tstamp - this->tstamp;
    std::cout << "Time lapse in Kalman: " << dt << "\n";
    _estimate(measurement, dt);
}

void KiteEKF::_estimate(const DM &measurement, const double &_dt)
{
    this->propagate(_dt);

    /** update step */
    DM y = measurement - DM::mtimes(H, State_Est);
    /** innovation covariance */
    DM S = DM::mtimes(DM::mtimes(H, Cov_Est), H.T()) + V;
    /** Kalman gain */
    DM K = DM::mtimes(DM::mtimes(Cov_Est, H.T()), DM::solve(S, SX::eye(7)));

    /** updated state estimate */
    DM x_est = State_Est + DM::mtimes(K, y);
    /** estimate covariance update */
    DM P_est = DM::mtimes((DM::eye(13) - DM::mtimes(K, H)), Cov_Est);

    State_Est = x_est;
    Cov_Est = P_est;
}
