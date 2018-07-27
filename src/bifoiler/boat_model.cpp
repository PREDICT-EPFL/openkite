#include <iostream> // TODO: remove later
#include <math.h>

#include <casadi/casadi.hpp>

#include "boat_model.h"
#include "boat_properties.h"
#include "boat_math.h"

using namespace casadi;
using namespace boat_math;

namespace bifoiler {

// TODO: debug only, remove later
DM x0, u0;
static void print_expr(const char *name, SX &sym, SX &x, SX &u)
{
    std::cout << name << "\n";
    std::cout << sym << "\n";
    Function f = Function(name,{x, u},{sym});
    std::cout << f(SXVector{x0, u0}) << "\n";
}

// TODO: debug only, remove later
static void print_jacobian(const char *name, SX &sym, SX &x, SX &u)
{
    SX jac = SX::jacobian(sym, x);
    print_expr(name, jac, x, u);
}

void BoatDynamics::Hydrodynamics(const SX &state, const SX &control, const BoatProperties &prop, SX &Fhbrf, SX &Mhbrf, SX &aoa, SX &ssa)
{
    double rho = prop.env.rho_sh2o;
    double AR  = prop.foils.ARff;        // [-], Aspect Ratio (Flügelstreckung)
    double S   = prop.foils.areaff;      // [m^2], Frontfoil Wing area
    double c   = prop.foils.mac;         // [m], Mean aerodynamic chord
    double b   = prop.foils.wingspanff;

    // Static Hydrodynamic Coefficient
    // All characteristics assumed linear
    double CL0         = prop.hydrodynamic.CL0;
    double CLa_tot     = prop.hydrodynamic.CLa_total;
    double e_o         = prop.hydrodynamic.e_oswald;
    double dw          = CLa_tot / (M_PI * e_o * AR); // downwash acting at the tail []
    double CD0_tot     = prop.hydrodynamic.CD0_total;
    double CYb         = prop.hydrodynamic.CYb;
    double Cm0         = prop.hydrodynamic.Cm0;
    double Cma         = prop.hydrodynamic.Cma;
    double Cn0         = prop.hydrodynamic.Cn0;
    double Cnb         = prop.hydrodynamic.Cnb;
    double Cl0         = prop.hydrodynamic.Cl0;
    double Clb         = prop.hydrodynamic.Clb;
    double CLq         = prop.hydrodynamic.CLq;
    double Cmq         = prop.hydrodynamic.Cmq;
    double CYr         = prop.hydrodynamic.CYr;
    double Cnr         = prop.hydrodynamic.Cnr;
    double Clr         = prop.hydrodynamic.Clr;
    double CYp         = prop.hydrodynamic.CYp;
    double Clp         = prop.hydrodynamic.Clp;
    double Cnp         = prop.hydrodynamic.Cnp;

    // Control Derivatives
    double CXdf        = prop.hydrodynamic.CXdf;
    double CYdr        = prop.hydrodynamic.CYdr;
    double CZde        = 0; // prop.hydrodynamic.CZde;
    double CZdf        = prop.hydrodynamic.CZdf;
    double Clda        = prop.hydrodynamic.CLda;
    double Cldr        = prop.hydrodynamic.CLdr;
    double Cmde        = 0; // prop.hydrodynamic.CMde;
    double Cmdf        = prop.hydrodynamic.CMdf;
    double Cnda        = prop.hydrodynamic.CNda;
    double Cndr        = prop.hydrodynamic.CNdr;

    // State variables definition
    SX v = state(Slice(0,3)); // Translational Velocity/Attitude towards Velocity Reference Frame
    SX W = state(Slice(3,6)); // Rotational Rate/Attitude Derivative

    // Control variables definition
    SX dF = control(0); // Flaps
    SX dA = control(1); // Aileron
    SX dR = control(2); // Rudder
    // SX dE = control(3); // Elevator

    // Variables and Orientation
    const double v1_eps   = 0.0001;    // Enable v(1)_0 = 0 w/o getting NaN
    SX V = SX::norm_2(v)+ v1_eps;        // Absolute Value of Velocity

    ssa = asin(v(1) / V);               // side slip angle [rad]
    aoa = atan2(v(2) , v(0) + v1_eps);  // angle of attack definition [rad]

    SX dyn_press = 0.5 * rho * pow(V,2);  // dynamic pressure

    SX q_aoa = SX::vertcat({cos(aoa/2), 0, sin(aoa/2), 0});
    SX q_ssa = SX::vertcat({cos(ssa/2), 0, 0, sin(-ssa/2)});
    SX q_BV = quatmul(q_aoa, q_ssa);

    // total drag coefficient
    SX CD = CD0_tot + pow(CL0 + CLa_tot * aoa,2) / (M_PI * e_o * AR);

    // Hydrodynamic Forces
    // Masterthesis equations. Assume Q = rho/2*S
    // Z  = (CL(aoa) + CLq * W(2) + CLaoa_dot(aoa_dot) + CLde*dE + CLdf*dF)*Q*V^2
    SX LIFT = (CL0 + CLa_tot * aoa + CZdf * dF) * dyn_press * S + (0.25 * CLq * c * S * rho) * V * W(1);

    // X  = (CD(aoa) + CXdf * dF) * Q * V^2
    SX DRAG = (CD + CXdf * dF) * dyn_press * S;

    // Y  = (CYb * ssa + CYp*W(1) + CYr* W(3) + CYdr * dR)*Q*V^2
    SX SF = (CYb * ssa + CYdr*dR) * dyn_press * S + 0.25 * (CYr * W(2) + CYp * W(0)) * (b * rho * S) * V;

    // Force return value
    Fhbrf = quatrot(q_BV, SX::vertcat({-DRAG, 0, -LIFT})) + SX::vertcat({0, SF, 0}); // why is SF already in brf?

    // Moments
    // Rolling Aerodynamic Moment
    // L = (Clb * ssa + Clp * W(1) + Clr * W(3) + Clda * dA + Cldr * dR) * Q * V^2 * b
    SX L = (Cl0 + Clb * ssa + Clda*dA + Cldr*dR) * dyn_press * S * b
           + (Clr * W(2) + Clp * W(0)) * (0.25 * rho * pow(b,2) * S) * V;

    // Pitching Aerodynamic Moment
    // M = (Cm0 + Cma*aoa + Cmaoa_dot*aoa_dot + Cmq*W(2) + Cmde*dE + Cmdf*dF) * Q * V^2* c
    SX M = (Cm0 + Cma * aoa  + Cmdf * dF) * dyn_press * S * c
           + Cmq * (0.25 * S * pow(c,2) * rho) * W(1) * V;

    // Yawing Aerodynamic Moment
    // N = (Cnb * ssa + Cnp * W(1) + Cnr * W(3) + Cnda * dA + Cndr * dR) * Q * V^2 * b
    SX N = (Cn0 + Cnb * ssa  + Cnda * dA + Cndr*dR) * dyn_press * S * b
           + (Cnp * W(0) + Cnr * W(2)) * (0.25 * S * pow(b,2) * rho) * V;

    // Angular motion equation in BRF
    SX Mhvrf = SX::vertcat({L, M, N});

    // Moment return value
    Mhbrf = quatrot(q_BV, Mhvrf);
}

void BoatDynamics::Hydrodynamics(const SX &state, const SX &control, const SX &params, const BoatProperties &prop, SX &Fhbrf, SX &Mhbrf, SX &aoa, SX &ssa)
{
    double rho = prop.env.rho_sh2o;
    double AR  = prop.foils.ARff;        // [-], Aspect Ratio (Flügelstreckung)
    double S   = prop.foils.areaff;      // [m^2], Frontfoil Wing area
    double c   = prop.foils.mac;         // [m], Mean aerodynamic chord
    double b   = prop.foils.wingspanff;

    SX CL0     = params(0);
    SX CLa_tot = params(1);
    double e_o         = prop.hydrodynamic.e_oswald;
    SX CD0_tot = params(2); //     = prop.hydrodynamic.CD0_total;
    SX CYb     = params(3); //    = prop.hydrodynamic.CYb;
    SX Cm0     = params(4); //    = prop.hydrodynamic.Cm0;
    SX Cma     = params(5); //    = prop.hydrodynamic.Cma;
    double Cn0 = prop.hydrodynamic.Cn0;
    SX Cnb     = params(6); //    = prop.hydrodynamic.Cnb;
    double Cl0 = prop.hydrodynamic.Cl0;
    SX Clb     = params(7); //    = prop.hydrodynamic.Clb;
    SX CLq     = params(8); //    = prop.hydrodynamic.CLq;
    SX Cmq     = params(9); //    = prop.hydrodynamic.Cmq;
    SX CYr     = params(10); //    = prop.hydrodynamic.CYr;
    SX Cnr     = params(11); //    = prop.hydrodynamic.Cnr;
    SX Clr     = params(12); //     = prop.hydrodynamic.Clr;
    SX CYp     = params(13); //    = prop.hydrodynamic.CYp;
    SX Clp     = params(14); //    = prop.hydrodynamic.Clp;
    SX Cnp     = params(15); //    = prop.hydrodynamic.Cnp;

    // Control Derivatives
    SX CXdf    = params(16); //       = prop.hydrodynamic.CXdf;
    SX CYdr    = params(17); //    = prop.hydrodynamic.CYdr;
    SX CZdf    = params(18); //    = prop.hydrodynamic.CZdf;
    SX Clda    = params(19); //    = prop.hydrodynamic.CLda;
    SX Cldr    = params(20); //    = prop.hydrodynamic.CLdr;
    SX Cmdf    = params(21); //    = prop.hydrodynamic.CMdf;
    SX Cnda    = params(22); //    = prop.hydrodynamic.CNda;
    SX Cndr    = params(23); //    = prop.hydrodynamic.CNdr;

    // State variables definition
    SX v = state(Slice(0,3)); // Translational Velocity/Attitude towards Velocity Reference Frame
    SX W = state(Slice(3,6)); // Rotational Rate/Attitude Derivative

    // Control variables definition
    SX dF = control(0); // Flaps
    SX dA = control(1); // Aileron
    SX dR = control(2); // Rudder
    // SX dE = control(3); // Elevator

    // Variables and Orientation
    const double v1_eps   = 0.0001;    // Enable v(1)_0 = 0 w/o getting NaN
    SX V = SX::norm_2(v)+ v1_eps;        // Absolute Value of Velocity

    ssa = asin(v(1) / V);               // side slip angle [rad]
    aoa = atan2(v(2) , v(0) + v1_eps);  // angle of attack definition [rad]

    SX dyn_press = 0.5 * rho * pow(V,2);  // dynamic pressure

    SX q_aoa = SX::vertcat({cos(aoa/2), 0, sin(aoa/2), 0});
    SX q_ssa = SX::vertcat({cos(ssa/2), 0, 0, sin(-ssa/2)});
    SX q_BV = quatmul(q_aoa, q_ssa);

    // total drag coefficient
    SX CD = CD0_tot + pow(CL0 + CLa_tot * aoa,2) / (M_PI * e_o * AR);

    // Hydrodynamic Forces
    // Masterthesis equations. Assume Q = rho/2*S
    // Z  = (CL(aoa) + CLq * W(2) + CLaoa_dot(aoa_dot) + CLde*dE + CLdf*dF)*Q*V^2
    SX LIFT = (CL0 + CLa_tot * aoa + CZdf * dF) * dyn_press * S + (0.25 * CLq * c * S * rho) * V * W(1);

    // X  = (CD(aoa) + CXdf * dF) * Q * V^2
    SX DRAG = (CD + CXdf * dF) * dyn_press * S;

    // Y  = (CYb * ssa + CYp*W(1) + CYr* W(3) + CYdr * dR)*Q*V^2
    SX SF = (CYb * ssa + CYdr*dR) * dyn_press * S + 0.25 * (CYr * W(2) + CYp * W(0)) * (b * rho * S) * V;

    // Force return value
    Fhbrf = quatrot(q_BV, SX::vertcat({-DRAG, 0, -LIFT})) + SX::vertcat({0, SF, 0}); // why is SF already in brf?

    // Moments
    // Rolling Aerodynamic Moment
    // L = (Clb * ssa + Clp * W(1) + Clr * W(3) + Clda * dA + Cldr * dR) * Q * V^2 * b
    SX L = (Cl0 + Clb * ssa + Clda*dA + Cldr*dR) * dyn_press * S * b
           + (Clr * W(2) + Clp * W(0)) * (0.25 * rho * pow(b,2) * S) * V;

    // Pitching Aerodynamic Moment
    // M = (Cm0 + Cma*aoa + Cmaoa_dot*aoa_dot + Cmq*W(2) + Cmde*dE + Cmdf*dF) * Q * V^2* c
    SX M = (Cm0 + Cma * aoa  + Cmdf * dF) * dyn_press * S * c
           + Cmq * (0.25 * S * pow(c,2) * rho) * W(1) * V;

    // Yawing Aerodynamic Moment
    // N = (Cnb * ssa + Cnp * W(1) + Cnr * W(3) + Cnda * dA + Cndr * dR) * Q * V^2 * b
    SX N = (Cn0 + Cnb * ssa  + Cnda * dA + Cndr*dR) * dyn_press * S * b
           + (Cnp * W(0) + Cnr * W(2)) * (0.25 * S * pow(b,2) * rho) * V;

    // Angular motion equation in BRF
    SX Mhvrf = SX::vertcat({L, M, N});

    // Moment return value
    Mhbrf = quatrot(q_BV, Mhvrf);
}

void BoatDynamics::Propulsion(const SX &state, const SX &control, const BoatProperties &prop, SX &Ftbrf, SX &Mtbrf)
{
    // TODO: move constants into BoatProperties

    SX thrust = control(3);

    // propmode: "endurance"
    // double Ftbrfx = 71;   // [N] in the Steady State
    // double Mtbrfx = 6.18; // [Nm], optimal torque (qopt). Prop rotates cw wrt brf_x => reaction

    // propmode: "speed"
    double Ftbrfx = 93;   // [N] in the Steady State
    double Mtbrfx = 4.77; // [Nm], optimal torque (qopt). Prop rotates cw wrt brf_x => reaction
    Ftbrf = SX::vertcat({Ftbrfx, 0, 0});
    Ftbrf = thrust * Ftbrf;

    // Parasite/Reaction Torque: Moment around y due to thrust
    // d_prop_vec = [0.6766 0.7716 0.8666]
    // boat_sims.cbdepth = 1
    // double d_prop_cog = d_prop_vec(boat_sims.cbdepth);
    double d_prop_cog = 0.7716;
    double Mrbrf = d_prop_cog*Ftbrfx;

    Mtbrf = SX::vertcat({Mtbrfx, Mrbrf, 0});
    Mtbrf = thrust * Mtbrf;
}


BoatDynamics::BoatDynamics(const BoatProperties &prop)
{
    double g = prop.env.g;
    double mass = prop.inertia.mass;
    double mass_cad = prop.inertia.mass_cad;
    double mcorr = mass / mass_cad; // Until better values: make boat uniformly heavier/inert

    double Ixy = mcorr * prop.inertia.Ixy;
    double Ixz = mcorr * prop.inertia.Ixz;
    double Iyz = mcorr * prop.inertia.Iyz;

    double Ixx = mcorr * prop.inertia.Ixx;
    double Iyy = mcorr * prop.inertia.Iyy;
    double Izz = mcorr * prop.inertia.Izz;

    SX Jbrf = SX::vertcat({  // [kg m^2], inertial tensor
        SX::horzcat({Ixx, Ixy, Ixz}),
        SX::horzcat({Ixy, Iyy, Iyz}),
        SX::horzcat({Ixz, Iyz, Izz})
    });

    // State and Control variables
    SX v    = SX::sym("v", 3);   // [m/s] translational velocity of the CoG in BRF
    SX W    = SX::sym("W", 3);   // [rad/s] angular velocities of boat in BRF
    SX r    = SX::sym("r", 3);   // [m] position of the CoG in IRF
    SX q_BI = SX::sym("q", 4);   // unit quaternion for transformation from IRF to BRF

    SX dF   = SX::sym("dF");     // Flaps
    SX dA   = SX::sym("dA");     // Aileron deflection [reserved, but not used] [rad]
    SX dR   = SX::sym("dR");     // Rudder deflection [rad]
    // SX dE   = SX::sym("dE");     // Elevator deflection [positive down] [rad]
    SX thrust = SX::sym("thrust"); // Thrust [-] (between 0.0 and 1.0)

    SX state    = SX::vertcat({v, W, r, q_BI});
    SX control  = SX::vertcat({dF, dA, dR, thrust}); // TODO: why is elevator deflection dE not used?

    // Gravity
    SX Fgbrf = mass*quatrot(q_BI, SX::vertcat({0,0,g})); // inertial reference frame, NED

    // Propulsion
    SX Ftbrf, Mtbrf;
    Propulsion(state, control, prop, Ftbrf, Mtbrf);

    // Hydrodynamics
    SX Fhbrf, Mhbrf, aoa, ssa; // TODO: aoa, ssa unused?
    Hydrodynamics(state, control, prop, Fhbrf, Mhbrf, aoa, ssa);

    // Buoyancy
    // SX Fbbrf, Mbbrf;
    // TODO

    // Damping
    // SX Df = SX::diag(SX::vertcat({-10, -20, -5})); // TODO: parametrize
    // SX Fdbrf = SX::mtimes(Df, v);

    // SX Dm = SX::diag(SX::vertcat({-2000, -800, -300}));
    // SX Mdbrf = SX::mtimes(Dm, W);


    SX Fbrf = Fhbrf + Ftbrf + Fgbrf; // + Fbbrf
    SX Mbrf = Mhbrf; // + Mtbrf + Mbbrf; // TODO: why no thrust moment?

    // Boat translational velocity in BRF
    SX v_dot_brf = Fbrf / mass - SX::cross(W,v);

    // Boat roataional velocity in BRF
    SX W_dot_brf = SX::mtimes(SX::inv(Jbrf), (Mbrf - SX::cross(W, SX::mtimes(Jbrf, W))));

    // Dynamic Equations: Kinematics
    SX r_dot_irf = quatrot_inverse(q_BI, v);
    SX q_BI_dot  = 0.5 * quatmul(q_BI, SX::vertcat({0,W})); // TODO: verify correctness

    // Differential equation
    SX dynamics = SX::vertcat({v_dot_brf, W_dot_brf, r_dot_irf, q_BI_dot});

    Function dynamics_func = Function("dynamics", {state, control}, {densify(dynamics)});
    SX jacobian = SX::jacobian(dynamics, state);
    Function jacobian_func = Function("jacobian", {state, control}, {densify(jacobian)});

    // define RK4 integrator scheme
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control.size1());
    SX dT = SX::sym("dT");

    // get symbolic expression for RK4 integrator
    SX integrator = rk4_symbolic(X, U, dynamics_func, dT);
    Function integrator_func = Function("integrator", {X,U,dT},{densify(integrator)});

    // assign class attributes
    this->State = state;
    this->Control = control;

    this->SymDynamics = dynamics;
    this->SymIntegartor = integrator;
    this->SymJacobian = jacobian;

    this->NumDynamics = dynamics_func;
    this->NumJacobian = jacobian_func;
    this->NumIntegrator = integrator_func;
}


/** identification constructor */
BoatDynamics::BoatDynamics(const BoatProperties &prop, const bool id_flag)
{
    double g = prop.env.g;
    double mass = prop.inertia.mass;
    double mass_cad = prop.inertia.mass_cad;
    double mcorr = mass / mass_cad; // Until better values: make boat uniformly heavier/inert

    double Ixy = mcorr * prop.inertia.Ixy;
    double Ixz = mcorr * prop.inertia.Ixz;
    double Iyz = mcorr * prop.inertia.Iyz;

    double Ixx = mcorr * prop.inertia.Ixx;
    double Iyy = mcorr * prop.inertia.Iyy;
    double Izz = mcorr * prop.inertia.Izz;

    SX Jbrf = SX::vertcat({  // [kg m^2], inertial tensor
        SX::horzcat({Ixx, Ixy, Ixz}),
        SX::horzcat({Ixy, Iyy, Iyz}),
        SX::horzcat({Ixz, Iyz, Izz})
    });

    // State and Control variables
    SX v    = SX::sym("v", 3);   // [m/s] translational velocity of the CoG in BRF
    SX W    = SX::sym("W", 3);   // [rad/s] angular velocities of boat in BRF
    SX r    = SX::sym("r", 3);   // [m] position of the CoG in IRF
    SX q_BI = SX::sym("q", 4);   // unit quaternion for transformation from IRF to BRF

    SX dF   = SX::sym("dF");     // Flaps
    SX dA   = SX::sym("dA");     // Aileron deflection [reserved, but not used] [rad]
    SX dR   = SX::sym("dR");     // Rudder deflection [rad]
    // SX dE   = SX::sym("dE");     // Elevator deflection [positive down] [rad]
    SX thrust = SX::sym("thrust"); // Thrust [-] (between 0.0 and 1.0)

    // model parameters
    SX CL0     = SX::sym("CL0");
    SX CLa_tot = SX::sym("CLa_tot");
    SX CD0_tot = SX::sym("CD0_tot"); //     = prop.hydrodynamic.CD0_total;
    SX CYb     = SX::sym("CYb"); //    = prop.hydrodynamic.CYb;
    SX Cm0     = SX::sym("Cm0"); //    = prop.hydrodynamic.Cm0;
    SX Cma     = SX::sym("Cma"); //    = prop.hydrodynamic.Cma;
    SX Cn0     = prop.hydrodynamic.Cn0;
    SX Cnb     = SX::sym("Cnb"); //    = prop.hydrodynamic.Cnb;
    SX Cl0     = prop.hydrodynamic.Cl0;
    SX Clb     = SX::sym("CLb"); //    = prop.hydrodynamic.Clb;
    SX CLq     = SX::sym("CLq"); //    = prop.hydrodynamic.CLq;
    SX Cmq     = SX::sym("Cmq"); //    = prop.hydrodynamic.Cmq;
    SX CYr     = SX::sym("CYr"); //    = prop.hydrodynamic.CYr;
    SX Cnr     = SX::sym("Cnr"); //    = prop.hydrodynamic.Cnr;
    SX Clr     = SX::sym("Clr"); //     = prop.hydrodynamic.Clr;
    SX CYp     = SX::sym("CYp"); //    = prop.hydrodynamic.CYp;
    SX Clp     = SX::sym("Clp"); //    = prop.hydrodynamic.Clp;
    SX Cnp     = SX::sym("Cnp"); //    = prop.hydrodynamic.Cnp;

    // Control Derivatives
    SX CXdf    = SX::sym("CXdf"); //       = prop.hydrodynamic.CXdf;
    SX CYdr    = SX::sym("CYdr"); //    = prop.hydrodynamic.CYdr;
    SX CZdf    = SX::sym("CZdf"); //    = prop.hydrodynamic.CZdf;
    SX Clda    = SX::sym("Clda"); //    = prop.hydrodynamic.CLda;
    SX Cldr    = SX::sym("Cldr"); //    = prop.hydrodynamic.CLdr;
    SX Cmdf    = SX::sym("Cmdf"); //    = prop.hydrodynamic.CMdf;
    SX Cnda    = SX::sym("Cnda"); //    = prop.hydrodynamic.CNda;
    SX Cndr    = SX::sym("Cndr"); //    = prop.hydrodynamic.CNdr;

    SX state    = SX::vertcat({v, W, r, q_BI});
    SX control  = SX::vertcat({dF, dA, dR, thrust}); // TODO: why is elevator deflection dE not used?
    SX params   = SX::vertcat({CL0, CLa_tot, CD0_tot, CYb, Cm0, Cma, Cnb, Clb, CLq, Cmq, CYr, Cnr,
                                Clr, CYp, Clp, Cnp, CXdf, CYdr, CZdf, Clda, Cldr, Cmdf, Cnda, Cndr});

    // Gravity
    SX Fgbrf = mass*quatrot(q_BI, SX::vertcat({0,0,g})); // inertial reference frame, NED

    // Propulsion
    SX Ftbrf, Mtbrf;
    Propulsion(state, control, prop, Ftbrf, Mtbrf);

    // Hydrodynamics
    SX Fhbrf, Mhbrf, aoa, ssa; // TODO: aoa, ssa unused - for some reason defined inside!
    Hydrodynamics(state, control, params, prop, Fhbrf, Mhbrf, aoa, ssa);

    // Buoyancy
    // SX Fbbrf, Mbbrf;
    // TODO

    // Damping
    //SX Df = SX::diag(SX::vertcat({-10, -20, -5})); // TODO: parametrize
    //SX Fdbrf = SX::mtimes(Df, v);

    //SX Dm = SX::diag(SX::vertcat({-2000, -800, -300}));
    //SX Mdbrf = SX::mtimes(Dm, W);


    SX Fbrf = Fhbrf + Ftbrf + Fgbrf; // + Fbbrf
    SX Mbrf = Mhbrf; // + Mtbrf + Mbbrf; // TODO: why no thrust moment?

    // Boat translational velocity in BRF
    SX v_dot_brf = Fbrf / mass - SX::cross(W,v);

    // Boat roataional velocity in BRF
    SX W_dot_brf = SX::mtimes(SX::inv(Jbrf), (Mbrf - SX::cross(W, SX::mtimes(Jbrf, W))));

    // Dynamic Equations: Kinematics
    SX r_dot_irf = quatrot_inverse(q_BI, v);
    SX q_BI_dot  = 0.5 * quatmul(q_BI, SX::vertcat({0,W})); // TODO: verify correctness

    // Differential equation
    SX dynamics = SX::vertcat({v_dot_brf, W_dot_brf, r_dot_irf, q_BI_dot});

    Function dynamics_func = Function("dynamics", {state, control, params}, {densify(dynamics)});
    SX jacobian = SX::jacobian(dynamics, state);
    Function jacobian_func = Function("jacobian", {state, control, params}, {densify(jacobian)});

    // define RK4 integrator scheme
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control.size1());
    SX dT = SX::sym("dT");

    // get symbolic expression for RK4 integrator
    //SX integrator = rk4_symbolic(X, U, dynamics_func, dT);
    //Function integrator_func = Function("integrator", {X,U,dT},{densify(integrator)});

    // assign class attributes
    this->State = state;
    this->Control = control;
    this->Params = params;

    this->SymDynamics = dynamics;
    //this->SymIntegartor = integrator;
    this->SymJacobian = jacobian;

    this->NumDynamics = dynamics_func;
    this->NumJacobian = jacobian_func;
    //this->NumIntegrator = integrator_func;
}

} // namespace bifoiler
