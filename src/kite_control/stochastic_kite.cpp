#include "stochastic_kite.hpp"

using namespace casadi;

StochasticKite::StochasticKite(const StochasticKiteProperties &props)
{
    SX theta  = SX::sym("theta");
    SX phi    = SX::sym("phi");
    SX gamma  = SX::sym("gamma");
    state = SX::vertcat({theta, phi, gamma});

    SX u   = SX::sym("u");
    control = u;

    double a = props.a;
    double b = props.b;
    double v = props.v;

    /** Dynamic equations */
    Dynamics = SX::vertcat({a * v * cos(theta) * cos(phi) * cos(gamma),
                            a * v * cos(theta) * sin(gamma),
                            b * v * cos(theta) * cos(phi) * u});
    NumDynamics = Function("Dynamics", {state, control}, {Dynamics});

    /** define output mapping */
    SX H = SX::zeros(2,3);
    H(0,0) = 1; H(1,1) = 1;
    OutputMap = Function("Map",{state}, {SX::mtimes(H, state)});
}

StochasticKite::StochasticKite()
{
    SX theta  = SX::sym("theta");
    SX phi    = SX::sym("phi");
    SX gamma  = SX::sym("gamma");
    state = SX::vertcat({theta, phi, gamma});

    SX u   = SX::sym("u");
    control = u;

    double a = 1.0;
    double b = 1.0;
    double v = 5.0;

    /** Dynamic equations */
    Dynamics = SX::vertcat({a * v * cos(theta) * cos(phi) * cos(gamma),
                            a * v * cos(theta) * sin(gamma),
                            b * v * cos(theta) * cos(phi) * u});
    NumDynamics = Function("Dynamics", {state, control}, {Dynamics});

    /** define output mapping */
    SX H = SX::zeros(2,3);
    H(0,0) = 1; H(1,1) = 1;
    OutputMap = Function("Map",{state}, {SX::mtimes(H, state)});
}

