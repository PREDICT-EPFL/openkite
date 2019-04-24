#ifndef STOCHASTIC_KITE_HPP
#define STOCHASTIC_KITE_HPP


#include "casadi/casadi.hpp"
#include <chrono>
#include "sys/stat.h"
#include <fstream>


struct StochasticKiteProperties
{
    StochasticKiteProperties(const double &_a = 1.0, const double &_b = 1.0, const double &_v = 5.0)
        : a(_a), b(_b), v(_v){}
    ~StochasticKiteProperties(){}

    double a, b, v;
};

/** Relatively simple mobile model : tricycle */
class StochasticKite
{
public:
    StochasticKite(const StochasticKiteProperties &props);
    StochasticKite();
    ~StochasticKite(){}

    casadi::Function getDynamics(){return NumDynamics;}
    casadi::Function getOutputMapping(){return OutputMap;}
private:
    casadi::SX state;
    casadi::SX control;
    casadi::SX Dynamics;

    casadi::Function NumDynamics;
    casadi::Function OutputMap;
};


#endif // STOCHASTIC_KITE_HPP
