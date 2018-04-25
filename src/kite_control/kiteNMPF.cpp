
#include "kiteNMPF.h"
#include "utility"
#include "pseudospectral/chebyshev.hpp"

using namespace casadi;

const DM KiteNMPF::DEFAULT_LBG = -DM::inf(1);
const DM KiteNMPF::DEFAULT_UBG = DM::inf(1);

const DM KiteNMPF::DEFAULT_LBX = -DM::inf(15);
const DM KiteNMPF::DEFAULT_UBX = DM::inf(15);

const DM KiteNMPF::DEFAULT_LBU = -DM::inf(4);
const DM KiteNMPF::DEFAULT_UBU = DM::inf(4);


KiteNMPF::KiteNMPF(std::shared_ptr<KiteDynamics> _Kite, const Function &_Path) : Kite(std::move(_Kite))
{
    PathFunc = _Path;

    /** set up default values */
    LBX = DEFAULT_LBX;
    UBX = DEFAULT_UBX;

    LBU = DEFAULT_LBU;
    UBU = DEFAULT_UBU;

    LBG = DEFAULT_LBG;
    UBG = DEFAULT_UBG;

    Q  =  0.5 * SX::diag(SX({1e1, 1e1, 1e2}));
    R  =  SX::diag(SX({1e-4, 1e-1, 1e-1, 1e-4}));
    W  = 1e-1;

    Scale_X = DM::eye(15); invSX = DM::eye(15);
    Scale_U = DM::eye(4);  invSU = DM::eye(4);

    NUM_SHOOTING_INTERVALS = 9;

    /** @attention : need sensible reference velocity */
    DM vel_ref = 0.05;
    this->setReferenceVelocity(vel_ref);

    WARM_START  = false;
    _initialized = false;
}


void KiteNMPF::createNLP()
{
    /** get dynamics function and state Jacobian */
    SX dynamics = Kite->getSymbolicDynamics();
    SX X = Kite->getSymbolicState();
    SX U = Kite->getSymbolicControl();

    /** state and control dimensionality */
    int n = 15;
    int m = 4;

    /** Order of polynomial interpolation */
    int N = NUM_SHOOTING_INTERVALS;

    /** define augmented dynamics of path parameter */
    SX V = SX::sym("V", 2);
    SX Uv = SX::sym("Uv");
    SX Av = SX::zeros(2,2); Av(0,1) = 1;
    SX Bv = SX::zeros(2,1); Bv(1,0) = 1;

    /** parameter dynamics */
    SX p_dynamics = SX::mtimes(Av, V) + SX::mtimes(Bv, Uv);

    /** augmented system */
    SX aug_state = SX::vertcat({X, V});
    SX aug_control = SX::vertcat({U, Uv});
    SX aug_dynamics = SX::vertcat({dynamics, p_dynamics});

    scale = true;

    /** evaluate augmented dynamics */
    Function aug_dynamo = Function("AUG_DYNAMO", {aug_state, aug_control}, {aug_dynamics});
    DynamicsFunc = aug_dynamo;

    /** ----------------------------------------------------------------------------------*/
    const int num_segments = 1;
    const int poly_order   = 9;
    const int dimx         = 15;
    const int dimu         = 4;
    const int dimp         = 0;

    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> spectral;
    SX diff_constr;

    SX x = SX::sym("x", dimx);
    SX u = SX::sym("u", dimu);

    if(scale)
    {
        SX SODE = aug_dynamo(SXVector{SX::mtimes(invSX,x), SX::mtimes(invSU, u)})[0];
        SODE = SX::mtimes(Scale_X, SODE);
        Function FunSODE = Function("scaled_ode", {x, u}, {SODE});

        diff_constr = spectral.CollocateDynamics(FunSODE, 0, 1.0);
    }
    else
    {
        diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, 1.0);
    }

    diff_constr = 0.1 * diff_constr;
    //diff_constr = diff_constr(Slice(0, diff_constr.size1() - dimx));

    /** define an integral cost */
    SX lagrange, residual;
    if(scale)
    {
        SXVector tmp = PathFunc(SXVector{SX::mtimes(invSX(13,13), x[13])});
        SX sym_path  = tmp[0];
        residual  = SX::mtimes(Scale_X(Slice(6,9), Slice(6,9)), sym_path) - x(Slice(6,9));
        lagrange  = SX::sumRows( SX::mtimes(Q, pow(residual, 2)) ) + SX::sumRows( SX::mtimes(W, pow(reference_velocity - x[14], 2)) );
        lagrange = lagrange + SX::sumRows( SX::mtimes(R, pow(u, 2)) );
    }
    else
    {
        SXVector tmp = PathFunc(SXVector{x[13]});
        SX sym_path  = tmp[0];
        residual  = sym_path - x(Slice(6,9));
        lagrange  = SX::sumRows( SX::mtimes(Q, pow(residual, 2)) ) + SX::sumRows( SX::mtimes(W, pow(reference_velocity - x[14], 2)) );
        lagrange = lagrange + SX::sumRows( SX::mtimes(R, pow(u, 2)) );
    }

    Function LagrangeTerm = Function("Lagrange", {x, u}, {lagrange});

    /** trace functions */
    PathError = Function("PathError", {x}, {residual});
    VelError  = Function("VelError", {x}, {reference_velocity - x[14]});

    SX mayer     =  SX::sumRows( SX::mtimes(Q, pow(residual, 2)) );
    Function MayerTerm    = Function("Mayer",{x}, {mayer});
    SX performance_idx = spectral.CollocateCost(MayerTerm, LagrangeTerm, 0, 1);

    SX varx = spectral.VarX();
    SX varu = spectral.VarU();

    SX opt_var = SX::vertcat(SXVector{varx, varu});

    /** debugging output */
    DynamicConstraints = Function("constraint_func", {opt_var}, {diff_constr});
    PerformanceIndex   = Function("performance_idx", {opt_var}, {performance_idx});

    SX lbg = SX::zeros(diff_constr.size());
    SX ubg = SX::zeros(diff_constr.size());

    /** set inequality (box) constraints */
    /** state */
    SX lbx = SX::repmat(SX::mtimes(Scale_X, LBX), poly_order + 1, 1);
    SX ubx = SX::repmat(SX::mtimes(Scale_X, UBX), poly_order + 1, 1);

    /** control */
    //lbx = SX::vertcat( {lbx, SX::repmat(SX::mtimes(Scale_U, LBU), poly_order, 1)} );
    //ubx = SX::vertcat( {ubx, SX::repmat(SX::mtimes(Scale_U, UBU), poly_order, 1)} );

    lbx = SX::vertcat( {lbx, SX::repmat(SX::mtimes(Scale_U, LBU), poly_order + 1, 1)} );
    ubx = SX::vertcat( {ubx, SX::repmat(SX::mtimes(Scale_U, UBU), poly_order + 1, 1)} );

    SX diff_constr_jacobian = SX::jacobian(diff_constr, opt_var);
    /** Augmented Jacobian */
    AugJacobian = Function("aug_jacobian",{opt_var}, {diff_constr_jacobian});

    /** formulate NLP */
    NLP["x"] = opt_var;
    NLP["f"] = performance_idx;
    NLP["g"] = diff_constr;

    OPTS["ipopt.linear_solver"]         = "ma97";
    OPTS["ipopt.print_level"]           = 0;
    OPTS["ipopt.tol"]                   = 1e-4;
    OPTS["ipopt.acceptable_tol"]        = 1e-4;
    OPTS["ipopt.max_iter"]              = 40;
    OPTS["ipopt.warm_start_init_point"] = "yes";
    NLP_Solver = nlpsol("solver", "ipopt", NLP, OPTS);

    /** set default args */
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    DM feasible_state = DM::mtimes(Scale_X, (UBX + LBX) / 2);
    DM feasible_control = DM::mtimes(Scale_U, (UBU + LBU) / 2);

    ARG["x0"] = DM::vertcat(DMVector{DM::repmat(feasible_state, poly_order + 1, 1),
                                     DM::repmat(feasible_control, poly_order + 1, 1)});
                                     //DM::repmat(feasible_control, poly_order, 1)});

    /** ----------------------------------------------------------------------------------*/

    /** NLP formulation */
    /** define optimisation variables */
    /** SX z = SX::sym("z", n, N+1);
    SX z_u = SX::sym("z_u", m, N); */

    /** allocate vectors to store constraints */
    /** SX vecx = {};
    SX g = {}; */

    /** state and control constraints */
    /** SX lbx = {};
    SX ubx = {}; */

    /** nonlinear state constraints */
    /** SX lbg = {};
    SX ubg = {}; */

    /** objective function */
    /** SX objective = 0; */

    /** Geometric path definition */
    /** SX theta = SX::sym("theta"); */

    /**  NLP formulation */
    /** Lagrangian term */
    /** SX pos = z(Slice(6,9), Slice(1, z.size2()-1));
    pos = SX::vec(pos); */

    /** path following error */
    /** DM ScaleXYZ = Scale_X(Slice(6,9), Slice(6,9));
    SX path_scaled = SX::mtimes(ScaleXYZ, SX::vertcat(PathFunc(SXVector{theta})));
    Function PathFuncScaled = Function("Scaled_Path",{theta},{path_scaled}); */

    /** @badcode: better DM to double conversion should exist */
    /** double scale_th = Scale_X(13,13).nonzeros()[0];
    SX path = kmath::mat_func( scale_th * z(13, Slice(1, z.size2()-1)), PathFuncScaled);
    SX err = pos - path;
    SX Qn = SX::kron(SX::eye(N-1), Q); */

    /** reference velocity following */
    /** SX err_v = z(14, Slice(1, z.size2()-1)).T() - SX::repmat(reference_velocity, N-1, 1);
    SX Qw = W * SX::eye(N-1);

    SX L = SX::sumRows( SX::mtimes(Qn, pow(err, 2)) ) + SX::sumRows( SX::mtimes(Qw, pow(err_v, 2)) ); */

    /** Mayer term */
    /** path = SX::vertcat(PathFuncScaled(SXVector{scale_th * z(13,0)}));
    SX err_n = z(Slice(6, 9), 0) - path;
    err_n = SX::vertcat( SXVector{err_n, z(14, 0) - reference_velocity});

    SX T = SX::mtimes((10 * SX::diagcat( {Q, W}) ), pow(err_n, 2)); */

    /** objective function */
    /** objective = L + SX::sumRows(T);

    CostFunction = Function("Cost", {z, z_u}, {objective});

    /** set equality constraints */
    /** differentiation matrix */
    /** DM _x, _D;
    std::pair<double, double> interval = std::make_pair<double, double>(0,1);
    kmath::cheb(_x, _D, N, interval);
    SX D = _D;
    D(D.size1()-1, Slice(0, D.size2())) = SX::zeros(N+1);
    D(D.size1() - 1, D.size2() - 1) = 1;
    SX Dn = SX::kron(D, SX::eye(n));

    SX iSx = SX::kron(SX::eye(N+1), invSX);
    SX iSu = SX::kron(SX::eye(N), invSU);

    SX F = kmath::mat_dynamics(SX::mtimes(invSX, z), SX::mtimes(invSU, z_u), aug_dynamo);
    SX G = SX::mtimes(SX::mtimes(Dn, iSx), SX::vec(z)) - F;
    G = G(Slice(0, N * n), 0);

    DynamicConstraints = Function("lox",{z, z_u},{G});

    lbg = SX::zeros(( N ) * n, 1);
    ubg = SX::zeros(( N ) * n, 1);

    /** set inequality (box) constraints */
    /** state */
    /** @todo: SCALE CONSTRAINTS */
    /** lbx = SX::repmat(SX::mtimes(Scale_X, LBX), N + 1, 1);
    ubx = SX::repmat(SX::mtimes(Scale_X, UBX), N + 1, 1);

    /** control */
     /** lbx = SX::vertcat( {lbx, SX::repmat(SX::mtimes(Scale_U, LBU), N, 1)} );
    ubx = SX::vertcat( {ubx, SX::repmat(SX::mtimes(Scale_U, UBU), N, 1)} );

    /** set decision variable */
    /** vecx = SX::vertcat({SX::vec(z), SX::vec(z_u)});

    /** formulate NLP */
    /** NLP["x"] = vecx;
    NLP["f"] = objective;
    NLP["g"] = G;

    OPTS["ipopt.linear_solver"]  = "ma97";
    OPTS["ipopt.print_level"]    = 0;
    OPTS["ipopt.tol"]            = 1e-5;
    OPTS["ipopt.acceptable_tol"] = 1e-4;
    //OPTS["ipopt.max_iter"]       = 15;
    NLP_Solver = nlpsol("solver", "ipopt", NLP, OPTS);

    /** set default args */
    /** ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    DM feasible_state = DM::mtimes(Scale_X, (UBX + LBX) / 2);
    DM feasible_control = DM::mtimes(Scale_U, (UBU + LBU) / 2);

    ARG["x0"] = DM::vertcat(DMVector{DM::repmat(feasible_state, N + 1, 1),
                                     DM::repmat(feasible_control, N, 1)}); */
}

void KiteNMPF::computeControl(const DM &_X0)
{
    int N = NUM_SHOOTING_INTERVALS;
    /** @badcode : remove magic constants */
    /** number of states */
    int nx = 15;
    /** number of control inputs */
    int nu = 4;

    /** rectify virtual state */
    DM X0 = _X0;
    //std::cout << "theta : before :" << X0[13] << " ";
    bool rectify = false;
    if(X0[13].nonzeros()[0] > 2 * M_PI)
    {
        X0[13] -= 2 * M_PI;
        rectify = true;
    }
    else if (X0[13].nonzeros()[0] < -2 * M_PI)
    {
        X0[13] += 2 * M_PI;
        rectify = true;
    }

    /** scale input */
    X0 = DM::mtimes(Scale_X, X0);
    double critical_val = Scale_X(13,13).nonzeros()[0] * 2 * M_PI;
    double flexibility  = Scale_X(13,13).nonzeros()[0] * 0.78;

    int idx_theta;

    if(WARM_START)
    {
        int idx_in = N * nx;
        int idx_out = idx_in + nx;
        ARG["lbx"](Slice(idx_in, idx_out), 0) = X0;
        ARG["ubx"](Slice(idx_in, idx_out), 0) = X0;

        /** relax virtual state constraint */
        idx_theta = idx_out - 2;
        ARG["lbx"](idx_theta) = X0[13] - flexibility;
        ARG["ubx"](idx_theta) = X0[13] + flexibility;

        ARG["lbx"](idx_theta + 1) = X0[14] - flexibility;
        ARG["ubx"](idx_theta + 1) = X0[14] + flexibility;

        /** rectify initial guess */
        if(rectify)
        {
            for(int i = 0; i < (N + 1) * nx; i += nx)
            {
                int idx = i + 13;
                if(NLP_X[idx].nonzeros()[0] > critical_val)
                    NLP_X[idx] -= critical_val;
                else if (NLP_X[idx].nonzeros()[0] < -critical_val)
                    NLP_X[idx] += critical_val;
            }
        }
        ARG["x0"]     = NLP_X;
        ARG["lam_g0"] = NLP_LAM_G;
        ARG["lam_x0"] = NLP_LAM_X;
    }
    else
    {
        ARG["x0"](Slice(0, (N + 1) * nx), 0) = DM::repmat(X0, (N + 1), 1);
        int idx_in = N * nx;
        int idx_out = idx_in + nx;
        ARG["lbx"](Slice(idx_in, idx_out), 0) = X0;
        ARG["ubx"](Slice(idx_in, idx_out), 0) = X0;

        /** relax virtual state constraint */
        idx_theta = idx_out - 2;
        ARG["lbx"](idx_theta) = X0[13] - flexibility;
        ARG["ubx"](idx_theta) = X0[13] + flexibility;

        ARG["lbx"](idx_theta + 1) = X0[14] - flexibility;
        ARG["ubx"](idx_theta + 1) = X0[14] + flexibility;
    }

    //DMVector dbg_prf    = PerformanceIndex({ARG["x0"]});
    //DMVector dbg_constr = DynamicConstraints({ARG["x0"]});

    //std::cout << "Cost Function : " << dbg_prf[0] << " Constraint Function : " << DM::norm_inf(dbg_constr[0]) << "\n";
    //DM state = ARG["x0"](Slice(N * nx, N * nx + nx));
    //std::cout << "State: " << DM::mtimes(invSX, state) << "\n";

    /** store optimal solution */
    DMDict res = NLP_Solver(ARG);
    NLP_X     = res.at("x");
    NLP_LAM_X = res.at("lam_x");
    NLP_LAM_G = res.at("lam_g");

    DM opt_x = NLP_X(Slice(0, (N + 1) * nx ));
    //DM invSX = DM::solve(Scale_X, DM::eye(15));
    OptimalTrajectory = DM::mtimes(invSX, DM::reshape(opt_x, nx, N + 1));
    DM opt_u = NLP_X( Slice((N + 1) * nx, NLP_X.size1()) );
    //DM invSU = DM::solve(Scale_U, DM::eye(4));
    OptimalControl = DM::mtimes(invSU, DM::reshape(opt_u, nu, N + 1));

    //std::cout << "Chosen : " << NLP_X[idx_theta] << "\n";

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
double KiteNMPF::getPathError()
{
    double error = 0;
    if(!OptimalTrajectory.is_empty())
    {
        DM state = OptimalTrajectory(Slice(0, OptimalTrajectory.size1()), OptimalTrajectory.size2() - 1);
        state = DM::mtimes(Scale_X, state);
        DMVector tmp = PathError(DMVector{state});
        error = DM::norm_2( tmp[0] ).nonzeros()[0];
    }
    return error;
}

/** get velocity error */
double KiteNMPF::getVelocityError()
{
    double error = 0;
    if(!OptimalTrajectory.is_empty())
    {
        DM state = OptimalTrajectory(Slice(0, OptimalTrajectory.size1()), OptimalTrajectory.size2() - 1);
        state = DM::mtimes(Scale_X, state);
        DMVector tmp = VelError(DMVector{state});
        error = DM::norm_2( tmp[0] ).nonzeros()[0];
    }
    return error;
}

/** get virtual state */
double KiteNMPF::getVirtState()
{
    double virt_state = 0;
    if(!OptimalTrajectory.is_empty())
    {
        virt_state = OptimalTrajectory(13, OptimalTrajectory.size2() - 1).nonzeros()[0];
    }
    return virt_state;
}

/** compute intial guess for virtual state */
DM KiteNMPF::findClosestPointOnPath(const DM &position, const DM &init_guess)
{
    SX theta = SX::sym("theta");
    SX sym_residual = 0.5 * SX::norm_2(PathFunc(SXVector{theta})[0] - SX(position));
    SX sym_gradient = SX::gradient(sym_residual,theta);

    Function grad = Function("gradient", {theta}, {sym_gradient});

    double tol = 1e-2;
    int counter = 0;
    DM theta_i = init_guess;
    DMVector result  = grad(DMVector{theta_i});
    double step = 0.25;

    /** check for local minimum */
    if(fabs(result[0].nonzeros()[0]) < tol)
    {
        theta_i = {M_PI_2 + 0.1};
        result  = grad(DMVector{theta_i});
    }

    /** @badcode : shitty code */
    while (fabs(result[0].nonzeros()[0]) >= tol)
    {
        counter++;
        theta_i -= step * grad(theta_i);
        /** gradient update */
        result  = grad(DMVector{theta_i});
        if(counter > 10)
            break;
    }
    std::cout << theta_i << "\n";
    return theta_i;
}

