#ifndef BOAT_MODEL_H
#define BOAT_MODEL_H

#include <casadi/casadi.hpp>
#include "boat_properties.h"

namespace bifoiler {

class BoatDynamics
{
public:
    BoatDynamics(const BoatProperties &prop);
    BoatDynamics(const BoatProperties &boat_prop, const bool id_flag);

    casadi::SX getSymbolicState(){return this->State;}
    casadi::SX getSymbolicControl(){return this->Control;}
    casadi::SX getSymbolicParams(){return this->Params;}

    casadi::SX getSymbolicDynamics(){return this->SymDynamics;}
    casadi::SX getSymbolicIntegrator(){return this->SymIntegartor;}
    casadi::SX getSymbolicJacobian(){return this->SymJacobian;}

    casadi::Function getNumericDynamics(){return this->NumDynamics;}
    casadi::Function getNumericIntegrator(){return this->NumIntegrator;}
    casadi::Function getNumericJacobian(){return this->NumJacobian;}

    static void Hydrodynamics(const casadi::SX &state,
                              const casadi::SX &control,
                              const casadi::SX &params,
                              const BoatProperties &prop,
                              casadi::SX &Fhbrf,
                              casadi::SX &Mhbrf,
                              casadi::SX &aoa,
                              casadi::SX &ssa);

    static void Hydrodynamics(const casadi::SX &state,
                              const casadi::SX &control,
                              const BoatProperties &prop,
                              casadi::SX &Fhbrf,
                              casadi::SX &Mhbrf,
                              casadi::SX &aoa,
                              casadi::SX &ssa);

    static void Propulsion(const casadi::SX &state,
                           const casadi::SX &control,
                           const BoatProperties &prop,
                           casadi::SX &Ftbrf,
                           casadi::SX &Mtbrf);
private:
    casadi::SX State;
    casadi::SX Control;
    casadi::SX Params;
    casadi::SX SymDynamics;
    casadi::SX SymIntegartor;
    casadi::SX SymJacobian;

    casadi::Function NumDynamics;
    casadi::Function NumIntegrator;
    casadi::Function NumJacobian;
};

} // namespace bifoiler

#endif /* BOAT_MODEL_H */
