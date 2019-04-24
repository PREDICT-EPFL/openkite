#ifndef NMPC_HPP
#define NMPC_HPP

#include <memory>
//#include "polymath.h"
#include "pseudospectral/chebyshev.hpp"

#define POLYMPC_USE_CONSTRAINTS

namespace polympc {

std::set<std::string> available_options = {"spectral.number_segments", "spectral.poly_order", "spectral.tf"};

template<typename ParamType>
ParamType get_param(const std::string &key, const casadi::Dict dict, const ParamType &default_val)
{
    /** find if the key is a valid option */
    if((available_options.find(key) != available_options.end()) && (dict.find(key) != dict.end()))
        return dict.find(key)->second;
    else if((available_options.find(key) == available_options.end()) && (dict.find(key) != dict.end()))
    {
        std::cout << "Unknown  parameter: " << key << "\n"; // << "Available parameters: " << available_options << "\n";
        return default_val;
    }
    else
        return default_val;
}

template <typename System, int NX, int NU, int NumSegments = 2, int PolyOrder = 5>
class nmpc
{
public:
    nmpc(const casadi::DM &_reference, const double &tf = 1.0, const casadi::DMDict &mpc_options = casadi::DMDict(), const casadi::Dict &solver_options = casadi::Dict());
    ~nmpc(){}

    /** contsraints setters */
    void setLBX(const casadi::DM &_lbx)
    {
        ARG["lbx"](casadi::Slice(0, NX * (PolyOrder * NumSegments + 1 ))) =
                   casadi::SX::repmat(casadi::SX::mtimes(Scale_X, _lbx), PolyOrder * NumSegments + 1, 1);
    }

    void setUBX(const casadi::DM &_ubx)
    {
        ARG["ubx"](casadi::Slice(0, NX * (PolyOrder * NumSegments + 1 ))) =
                   casadi::SX::repmat(casadi::SX::mtimes(Scale_X, _ubx), PolyOrder * NumSegments + 1, 1);
    }

    void setLBU(const casadi::DM &_lbu)
    {
        int start = NX * (PolyOrder * NumSegments + 1 );
        int finish = start + NU * (PolyOrder * NumSegments + 1 );
        ARG["lbx"](casadi::Slice(start, finish)) = casadi::SX::repmat(casadi::SX::mtimes(Scale_U, _lbu), PolyOrder * NumSegments + 1, 1);
    }
    void setUBU(const casadi::DM &_ubu)
    {
        int start = NX * (PolyOrder * NumSegments + 1 );
        int finish = start + NU * (PolyOrder * NumSegments + 1 );
        ARG["ubx"](casadi::Slice(start, finish)) = casadi::SX::repmat(casadi::SX::mtimes(Scale_U, _ubu), PolyOrder * NumSegments + 1, 1);
    }

    void setStateScaling(const casadi::DM &Scaling){Scale_X = Scaling;
                                                      invSX = casadi::DM::solve(Scale_X, casadi::DM::eye(Scale_X.size1()));}
    void setControlScaling(const casadi::DM &Scaling){Scale_U = Scaling;
                                                      invSU = casadi::DM::solve(Scale_U, casadi::DM::eye(Scale_U.size1()));}

    void createNLP(const casadi::Dict &solver_options);
    void updateParams(const casadi::Dict &params);

    void enableWarmStart(){WARM_START = true;}
    void disableWarmStart(){WARM_START = false;}
    void computeControl(const casadi::DM &_X0);

    casadi::DM getOptimalControl(){return OptimalControl;}
    casadi::DM getOptimalTrajetory(){return OptimalTrajectory;}

    casadi::Dict getStats(){return stats;}
    bool initialized(){return _initialized;}

    double getPathError();

private:
    System system;
    casadi::SX Reference;
    uint   nx, nu, ny, np;
    double Tf;

    casadi::SX       Contraints;
    casadi::Function ContraintsFunc;

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
    casadi::SX Q, R, W;

    casadi::DM NLP_X, NLP_LAM_G, NLP_LAM_X;
    casadi::Function NLP_Solver;
    casadi::SXDict NLP;
    casadi::Dict OPTS;
    casadi::DMDict ARG;
    casadi::Dict stats;

    casadi::DM OptimalControl;
    casadi::DM OptimalTrajectory;

    unsigned NUM_COLLOCATION_POINTS;
    bool WARM_START;
    bool _initialized;
    bool scale;

    /** TRACE FUNCTIONS */
    casadi::Function DynamicsFunc;
    casadi::Function DynamicConstraints;
    casadi::Function PerformanceIndex;
    casadi::Function CostFunction;
    casadi::Function PathError;

    casadi::Function m_Jacobian;
    casadi::Function m_Dynamics;
};

template<typename System, int NX, int NU, int NumSegments, int PolyOrder>
nmpc<System, NX, NU, NumSegments, PolyOrder>::nmpc(const casadi::DM &_reference, const double &tf, const casadi::DMDict &mpc_options, const casadi::Dict &solver_options)
{
    /** set up default */
    casadi::Function dynamics = system.getDynamics();
    nx = dynamics.nnz_out();
    nu = dynamics.nnz_in() - nx;
    Tf = tf;

    assert(NX == nx);
    assert(NU == nu);

    casadi::Function output   = system.getOutputMapping();
    ny = output.nnz_out();
    Reference = _reference;

    assert(ny == Reference.size1());

    Q  = casadi::SX::diag(casadi::SX({10,150,50}));
    R  = casadi::SX::diag(casadi::SX({0.1,1,0.1}));

    /** problem scaling */
    scale = false;
    if(mpc_options.find("mpc.scaling") != mpc_options.end())
        scale = static_cast<bool>(mpc_options.find("mpc.scaling")->second.nonzeros()[0]);

    if(mpc_options.find("mpc.scale_x") != mpc_options.end())
    {
        Scale_X = mpc_options.find("mpc.scale_x")->second;
        assert(NX == Scale_X.size1());
        assert(NX == Scale_X.size2());
        invSX = casadi::DM::solve(Scale_X, casadi::DM::eye(Scale_X.size1()));
    }

    if(mpc_options.find("mpc.scale_u") != mpc_options.end())
    {
        Scale_U = mpc_options.find("mpc.scale_u")->second;
        assert(NU == Scale_U.size1());
        assert(NU == Scale_U.size2());
        invSU = casadi::DM::solve(Scale_U, casadi::DM::eye(Scale_U.size1()));
    }

    /** assume unconstrained problem */
    LBX = -casadi::DM::inf(nx);
    UBX = casadi::DM::inf(nx);
    LBU = -casadi::DM::inf(nu);
    UBU = casadi::DM::inf(nu);

    WARM_START  = false;
    _initialized = false;

    /** create NLP */
    createNLP(solver_options);
}

/** update solver paramters */
template<typename System, int NX, int NU, int NumSegments, int PolyOrder>
void nmpc<System, NX, NU, NumSegments, PolyOrder>::updateParams(const casadi::Dict &params)
{
    for (casadi::Dict::const_iterator it = params.begin(); it != params.end(); ++it)
    {
        OPTS[it->first] = it->second;
    }
}

template<typename System, int NX, int NU, int NumSegments, int PolyOrder>
void nmpc<System, NX, NU, NumSegments, PolyOrder>::createNLP(const casadi::Dict &solver_options)
{
    /** get dynamics function and state Jacobian */
    casadi::Function dynamics = system.getDynamics();
    casadi::Function output   = system.getOutputMapping();
    casadi::SX x = casadi::SX::sym("x", nx);
    casadi::SX u = casadi::SX::sym("u", nu);
    DynamicsFunc = dynamics;

    if(scale)
    {
        assert(Scale_X.size1() != 0);
        assert(Scale_U.size1() != 0);
    }

    /** ----------------------------------------------------------------------------------*/
    /** set default properties of approximation */
    const int num_segments = NumSegments;  //get_param<int>("spectral.number_segments", spectral_props.props, 2);
    const int poly_order   = PolyOrder;    //get_param<int>("spectral.poly_order", spectral_props.props, 5);
    const int dimx         = NX;
    const int dimu         = NU;
    const int dimp         = 0;
    const double tf        = Tf;           //get_param<double>("spectral.tf", spectral_props.props, 1.0);

    NUM_COLLOCATION_POINTS = num_segments * poly_order;
    /** Order of polynomial interpolation */

    Chebyshev<casadi::SX, poly_order, num_segments, dimx, dimu, dimp> spectral;
    casadi::SX diff_constr;

    if(scale)
    {
        casadi::SX SODE = dynamics(casadi::SXVector{casadi::SX::mtimes(invSX, x), casadi::SX::mtimes(invSU, u)})[0];
        SODE = casadi::SX::mtimes(Scale_X, SODE);
        casadi::Function FunSODE = casadi::Function("scaled_ode", {x, u}, {SODE});

        diff_constr = spectral.CollocateDynamics(FunSODE, 0, tf);
    }
    else
    {
        diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, tf);
    }

    diff_constr = diff_constr(casadi::Slice(0, diff_constr.size1() - dimx));

    /** define an integral cost */
    casadi::SX lagrange, residual;
    if(scale)
    {
        casadi::SX _invSX = invSX(casadi::Slice(0, NX), casadi::Slice(0, NX));
        residual  = Reference - output({casadi::SX::mtimes(_invSX, x)})[0];
        lagrange  = casadi::SX::sumRows( casadi::SX::mtimes(Q, pow(residual, 2)) );
        lagrange = lagrange + casadi::SX::sumRows( casadi::SX::mtimes(R, pow(u, 2)) );
    }
    else
    {
        residual  = Reference - output({x})[0];
        lagrange  = casadi::SX::sumRows( casadi::SX::mtimes(Q, pow(residual, 2)) );
        lagrange = lagrange + casadi::SX::sumRows( casadi::SX::mtimes(R, pow(u, 2)) );
    }

    casadi::Function LagrangeTerm = casadi::Function("Lagrange", {x, u}, {lagrange});

    /** trace functions */
    PathError = casadi::Function("PathError", {x}, {residual});

    casadi::SX mayer           =  casadi::SX::sumRows( casadi::SX::mtimes(Q, pow(residual, 2)) );
    casadi::Function MayerTerm = casadi::Function("Mayer",{x}, {mayer});
    casadi::SX performance_idx = spectral.CollocateCost(MayerTerm, LagrangeTerm, 0.0, tf);

    casadi::SX varx = spectral.VarX();
    casadi::SX varu = spectral.VarU();

    casadi::SX opt_var = casadi::SX::vertcat(casadi::SXVector{varx, varu});

    /** debugging output */
    DynamicConstraints = casadi::Function("constraint_func", {opt_var}, {diff_constr});
    PerformanceIndex   = casadi::Function("performance_idx", {opt_var}, {performance_idx});

    casadi::SX lbg = casadi::SX::zeros(diff_constr.size());
    casadi::SX ubg = casadi::SX::zeros(diff_constr.size());

    /** set inequality (box) constraints */
    /** state */
    casadi::SX lbx = casadi::SX::repmat(casadi::SX::mtimes(Scale_X, LBX), poly_order * num_segments + 1, 1);
    casadi::SX ubx = casadi::SX::repmat(casadi::SX::mtimes(Scale_X, UBX), poly_order * num_segments + 1, 1);

    /** control */
    lbx = casadi::SX::vertcat( {lbx, casadi::SX::repmat(casadi::SX::mtimes(Scale_U, LBU), poly_order * num_segments + 1, 1)} );
    ubx = casadi::SX::vertcat( {ubx, casadi::SX::repmat(casadi::SX::mtimes(Scale_U, UBU), poly_order * num_segments + 1, 1)} );

    /** experimental : constraints*/
#ifdef POLYMPC_USE_CONSTRAINTS

    /** bounding box vertices */
    casadi::SX p1 = casadi::SX::vertcat({-2.5,-1});
    casadi::SX p2 = casadi::SX::vertcat({-2.5, 1});
    casadi::SX p3 = casadi::SX::vertcat({-7,  -1});
    casadi::SX p4 = casadi::SX::vertcat({-7,  1});

    /** second bounding box */
    casadi::SX q1 = casadi::SX::vertcat({2.5,-1});
    casadi::SX q2 = casadi::SX::vertcat({2.5, 1});
    casadi::SX q3 = casadi::SX::vertcat({7,  -1});
    casadi::SX q4 = casadi::SX::vertcat({7,  1});

    /** robot vertices */
    casadi::SX r1 = casadi::SX::vertcat({1.2, -0.6});
    casadi::SX r2 = casadi::SX::vertcat({1.2,  0.6});
    casadi::SX r3 = casadi::SX::vertcat({-1.2,  0.6});
    casadi::SX r4 = casadi::SX::vertcat({-1.2, -0.6});

    /** define constraints function */
    casadi::SX p = x(casadi::Slice(3,5));
    casadi::Function th_f = casadi::Function("th_f",{x,u},{x[5]});
    casadi::Function xy_f = casadi::Function("xy_f",{x,u}, {p});
    casadi::SX xy = spectral.CollocateFunction(xy_f);
    casadi::SX th = spectral.CollocateFunction(th_f);

    int n_steps = PolyOrder * NumSegments + 1;
    casadi::SX n1 = casadi::SX::sym("n", 2 * n_steps); // plane1 direction vector
    casadi::SX b1 = casadi::SX::sym("b", n_steps);     // plane1 bias

    casadi::SX n2 = casadi::SX::sym("n", 2 * n_steps); // plane2 direction vector
    casadi::SX b2 = casadi::SX::sym("b", n_steps);     // plane2 bias
    casadi::SX ineq_constr;

    /** set up the hyperplane constraints */
    for(int i = 0; i < n_steps; ++i)
    {
        int idx = i * 2;
        casadi::Slice slice_i = casadi::Slice(idx, idx + 2);
        /** compute the rotation matrix */
        casadi::SX th_i = th[i];
        casadi::SX col1 = casadi::SX::vertcat({cos(th_i), sin(th_i)});
        casadi::SX col2 = casadi::SX::vertcat({-sin(th_i), cos(th_i)});
        casadi::SX RX = casadi::SX::horzcat({col1,col2});
        /** find root points wrt to parking place */
        casadi::SX r1_w = casadi::SX::mtimes(RX, r1) + xy(slice_i);
        casadi::SX r2_w = casadi::SX::mtimes(RX, r2) + xy(slice_i);
        casadi::SX r3_w = casadi::SX::mtimes(RX, r3) + xy(slice_i);
        casadi::SX r4_w = casadi::SX::mtimes(RX, r4) + xy(slice_i);


        /** search for the first separating hyperplane */
        casadi::SX cr1 = casadi::SX::dot(-n1(slice_i), r1_w) + b1[i];
        casadi::SX cr2 = casadi::SX::dot(-n1(slice_i), r2_w) + b1[i];
        casadi::SX cr3 = casadi::SX::dot(-n1(slice_i), r3_w) + b1[i];
        casadi::SX cr4 = casadi::SX::dot(-n1(slice_i), r4_w) + b1[i];

        casadi::SX c2 = casadi::SX::dot(n1(slice_i), p1) - b1[i];
        casadi::SX c3 = casadi::SX::dot(n1(slice_i), p2) - b1[i];
        casadi::SX c4 = casadi::SX::dot(n1(slice_i), p3) - b1[i];
        casadi::SX c5 = casadi::SX::dot(n1(slice_i), p4) - b1[i];
        casadi::SX c6 = casadi::SX::dot(n1(slice_i), n1(slice_i)) - 1;

        ineq_constr = casadi::SX::vertcat({ineq_constr, cr1, cr2, cr3, cr4, c2, c3, c4, c5, c6});

        /** search for the second separating hyperplane */
        casadi::SX cr1_2 = casadi::SX::dot(-n2(slice_i), r1_w) + b2[i];
        casadi::SX cr2_2 = casadi::SX::dot(-n2(slice_i), r2_w) + b2[i];
        casadi::SX cr3_2 = casadi::SX::dot(-n2(slice_i), r3_w) + b2[i];
        casadi::SX cr4_2 = casadi::SX::dot(-n2(slice_i), r4_w) + b2[i];

        casadi::SX c2_2 = casadi::SX::dot(n2(slice_i), q1) - b2[i];
        casadi::SX c3_2 = casadi::SX::dot(n2(slice_i), q2) - b2[i];
        casadi::SX c4_2 = casadi::SX::dot(n2(slice_i), q3) - b2[i];
        casadi::SX c5_2 = casadi::SX::dot(n2(slice_i), q4) - b2[i];
        casadi::SX c6_2 = casadi::SX::dot(n2(slice_i), n2(slice_i)) - 1;

        ineq_constr = casadi::SX::vertcat({ineq_constr, cr1_2, cr2_2, cr3_2, cr4_2, c2_2, c3_2, c4_2, c5_2, c6_2});
    }

    /** remove the initial point check */
    //ineq_constr = ineq_constr(casadi::Slice(0, ineq_constr.size1() - xy_f.nnz_out()));

    /** define the inequality bounds*/
    casadi::SX upper = casadi::SX::zeros(ineq_constr.size1());
    casadi::SX lower = -casadi::SX::inf(ineq_constr.size1());

    /** append G constraint */
    diff_constr = casadi::SX::vertcat(casadi::SXVector{diff_constr, ineq_constr});
    lbg = casadi::SX::vertcat({lbg, lower});
    ubg = casadi::SX::vertcat({ubg, upper});

    /** append optimization vector */
    opt_var = casadi::SX::vertcat({opt_var, n1, b1, n2, b2});
    casadi::SX plane_x0 = casadi::SX::zeros(2 * 3 * n_steps);
    casadi::SX plane_lb = -casadi::SX::inf(2 * 3 * n_steps);
    casadi::SX plane_ub = casadi::SX::inf(2 * 3 * n_steps);

    lbx = casadi::SX::vertcat({lbx, plane_lb});
    ubx = casadi::SX::vertcat({ubx, plane_ub});

    std::cout << "Use constraints \n";

#endif

    casadi::SX diff_constr_jacobian = casadi::SX::jacobian(diff_constr, opt_var);
    /** Augmented Jacobian */
    m_Jacobian = casadi::Function("aug_jacobian",{opt_var}, {diff_constr_jacobian});

    /** formulate NLP */
    NLP["x"] = opt_var;
    NLP["f"] = performance_idx; //  1e-3 * casadi::SX::dot(diff_constr, diff_constr);
    NLP["g"] = diff_constr;

    /** default solver options */
    OPTS["ipopt.linear_solver"]         = "ma97";
    OPTS["ipopt.print_level"]           = 1;
    OPTS["ipopt.tol"]                   = 1e-4;
    OPTS["ipopt.acceptable_tol"]        = 1e-4;
    OPTS["ipopt.max_iter"]              = 150;
    OPTS["ipopt.warm_start_init_point"] = "yes";
    //OPTS["ipopt.hessian_approximation"] = "limited-memory";

    /** set user defined options */
    if(!solver_options.empty())
        updateParams(solver_options);

    NLP_Solver = casadi::nlpsol("solver", "ipopt", NLP, OPTS);

    /** set default args */
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    casadi::DM feasible_state = casadi::DM::zeros(UBX.size());
    casadi::DM feasible_control = casadi::DM::zeros(UBU.size());

    ARG["x0"] = casadi::DM::vertcat(casadi::DMVector{casadi::DM::repmat(feasible_state, poly_order * num_segments + 1, 1),
                                     casadi::DM::repmat(feasible_control, poly_order * num_segments + 1, 1), plane_x0});
}

template<typename System, int NX, int NU, int NumSegments, int PolyOrder>
void nmpc<System, NX, NU, NumSegments, PolyOrder>::computeControl(const casadi::DM &_X0)
{
    int N = NUM_COLLOCATION_POINTS;

    /** rectify virtual state */
    casadi::DM X0 = casadi::DM::mtimes(Scale_X, _X0);

    /** scale input */
    int idx_theta;
    std::cout << "Compute control at: " << X0 << "\n";

    if(WARM_START)
    {
        int idx_in = N * NX;
        int idx_out = idx_in + NX;
        ARG["lbx"](casadi::Slice(idx_in, idx_out), 0) = X0;
        ARG["ubx"](casadi::Slice(idx_in, idx_out), 0) = X0;

        ARG["x0"]     = NLP_X;
        ARG["lam_g0"] = NLP_LAM_G;
        ARG["lam_x0"] = NLP_LAM_X;
    }
    else
    {
        ARG["x0"](casadi::Slice(0, (N + 1) * NX), 0) = casadi::DM::repmat(X0, (N + 1), 1);
        int idx_in = N * NX;
        int idx_out = idx_in + NX;
        ARG["lbx"](casadi::Slice(idx_in, idx_out), 0) = X0;
        ARG["ubx"](casadi::Slice(idx_in, idx_out), 0) = X0;
    }

    /** store optimal solution */
    casadi::DMDict res = NLP_Solver(ARG);
    NLP_X     = res.at("x");
    NLP_LAM_X = res.at("lam_x");
    NLP_LAM_G = res.at("lam_g");

    casadi::DM opt_x = NLP_X(casadi::Slice(0, (N + 1) * NX));
    //DM invSX = DM::solve(Scale_X, DM::eye(15));
    OptimalTrajectory = casadi::DM::mtimes(invSX, casadi::DM::reshape(opt_x, NX, N + 1));
    //casadi::DM opt_u = NLP_X( casadi::Slice((N + 1) * NX, NLP_X.size1()) );
    casadi::DM opt_u = NLP_X( casadi::Slice((N + 1) * NX, (N + 1) * NX + (N + 1) * NU ) );
    //DM invSU = DM::solve(Scale_U, DM::eye(4));
    OptimalControl = casadi::DM::mtimes(invSU, casadi::DM::reshape(opt_u, NU, N + 1));

    stats = NLP_Solver.stats();
    std::cout << stats << "\n";

    std::string solve_status = static_cast<std::string>(stats["return_status"]);
    if(solve_status.compare("Invalid_Number_Detected") == 0)
    {
        std::cout << "X0 : " << ARG["x0"] << "\n";
        //assert(false);
    }
    if(solve_status.compare("Infeasible_Problem_Detected") == 0)
    {
        std::cout << "X0 : " << ARG["x0"] << "\n";
        //assert(false);
    }

    enableWarmStart();
}

/** get path error */
template<typename System, int NX, int NU, int NumSegments, int PolyOrder>
double nmpc<System, NX, NU, NumSegments, PolyOrder>::getPathError()
{
    double error = 0;
    if(!OptimalTrajectory.is_empty())
    {
        casadi::DM state = OptimalTrajectory(casadi::Slice(0, OptimalTrajectory.size1()), OptimalTrajectory.size2() - 1);
        state = casadi::DM::mtimes(Scale_X, state);
        casadi::DMVector tmp = PathError(casadi::DMVector{state});
        error = casadi::DM::norm_2( tmp[0] ).nonzeros()[0];
    }
    return error;
}


} //polympc namespace


#endif // NMPC_HPP
