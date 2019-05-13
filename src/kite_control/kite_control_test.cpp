#include "kiteEKF.h"
#include "kiteNMPF.h"

#define BOOST_TEST_TOOLS_UNDER_DEBUGGER
#define BOOST_TEST_MODULE kite_control_test
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include "pseudospectral/chebyshev.hpp"
#include <unordered_set>

using namespace casadi;

BOOST_AUTO_TEST_SUITE( kite_control_suite_test )

BOOST_AUTO_TEST_CASE( integrator_test )
{
    AlgorithmProperties algo_props;
    algo_props.Integrator = CVODES;
    algo_props.sampling_time = 0.02;

    /** create a rigid body object */
    RigidBodyKinematics rigid_body = RigidBodyKinematics(algo_props);
    Function fIntegrator = rigid_body.getNumericIntegrator();
    Function fJacobian = rigid_body.getNumericJacobian();
    Function fDynamics = rigid_body.getNumericDynamcis();

    /** perform integration step */
    DM init_state = DM::vertcat({4.318732, 0.182552, 0.254833, 1.85435, -0.142882, -0.168359,
                              -0.229383, -0.0500282, -0.746832, 0.189409, -0.836349, -0.48178, 0.180367});

    DM dynamics = fDynamics(DMVector{init_state});
    std::cout << "Dynamics: " << dynamics(0) << "\n";

    DM control = DM::zeros(3);
    //DMDict args = {{"x0", init_state}, {"p", control}};
    DMDict args = {{"x0", init_state}};
    DMDict out = fIntegrator(args);
    DM new_state = out["xf"];

    std::cout << "Integration: init_state: " << init_state << "\n";
    std::cout << "Integration: new_state: " << new_state << "\n" << "\n";

    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE( ekf_test )
{
    /** make test data */
    double dt = 0.0084;
    DM control = DM::vertcat({0,0,0});
    DM measurement = DM::vertcat({1.4522, -3.1274, -1.7034, -0.5455, -0.2382, -0.2922, -0.7485});
    DM x_est = DM::vertcat({6.0026, -0.3965, 0.1705,
                        0.4414, -0.2068, 0.9293,
                        1.4634, -3.1765, -1.7037,
                       -0.5486, -0.2354, -0.2922, -0.7471});
    /** reference estimation from matlab */
    DM reference_est = DM::vertcat({5.9982, -0.3819, 0.1637,
                                    0.3578, -0.1900, 0.8774,
                                    1.4522, -3.1274, -1.7034,
                                    -0.5455, -0.2382, -0.2922, -0.7485});

    /** create a filter and perform estimation */
    std::string kite_config_file = "umx_radian.yaml";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_config_file);
    AlgorithmProperties algo_props;
    algo_props.Integrator = RK4;

    KiteEKF estimator = KiteEKF(kite_props, algo_props);
    estimator.setControl(control);
    estimator.setEstimation(x_est);

    std::chrono::time_point<std::chrono::system_clock> start = kite_utils::get_time();

    estimator._estimate(measurement, dt);
    x_est = estimator.getEstimation();
    //estimator.propagate(dt);

    std::chrono::time_point<std::chrono::system_clock> stop = kite_utils::get_time();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    auto latency = static_cast<double>(duration.count());

    std::cout << "EKF_TEST COMPUTATION TIME: " << latency * 1e-6 << " [seconds]" << "\n";

    //BOOST_CHECK(DM::norm_inf(x_est - reference_est).nonzeros()[0] < 0.01);
    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_CASE( lqr_test )
{
    /** lyapunov test */
    Eigen::MatrixXd A(2,2), Q(2,2);
    A << 1, 2, -3, -4;
    Q << 3, 1, 1, 1;
    kite_utils::time_point start = kite_utils::get_time();
    Eigen::MatrixXd X = kmath::oc::lyapunov(A, Q);
    kite_utils::time_point finish = kite_utils::get_time();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    auto latency = static_cast<double>(duration.count());

    std::cout << "Lyapunov equation solution: \n" << X << "\n";
    std::cout << (A * X + X * A.transpose()) << "\n";
    std::cout << "Time elapsed : " << latency * 1e-6 << " [seconds]" << "\n";

    /** PINV test */
    start = kite_utils::get_time();
    Eigen::MatrixXd pinvA = kmath::oc::pinv(A);
    finish = kite_utils::get_time();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    latency = static_cast<double>(duration.count());
    std::cout << "pinv: " << pinvA << "\n";
    std::cout << "Time elapsed : " << latency * 1e-6 << " [seconds]" << "\n";

    /** INIT Newton test */
    X = kmath::oc::init_newton_care(A, Q);
    std::cout << "Init: \n" << X << "\n";

    Eigen::MatrixXd A1(2,2), B1(2,2), C1(2,2);
    A1 << -3, 2, 1, 1;
    B1 << 0, 0, 0, 1;
    C1 << 1, -1, -1, 1;
    start = kite_utils::get_time();
    Eigen::MatrixXd S = kmath::oc::care(A1, B1, C1);
    finish = kite_utils::get_time();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    latency = static_cast<double>(duration.count());
    std::cout << "CARE solved: \n" << S << "\n";
    std::cout << "Time elapsed : " << latency * 1e-6 << " [seconds]" << "\n";

    std::cout <<"-----------------------Controllability test------------------------\n";
    Eigen::MatrixXd Ac(2,2), Bc(2,2), Qc(2,2), R(2,2), M(2,2);
    Ac << 1,1,4,-2;
    Bc << Eigen::MatrixXd::Identity(2,2);
    Qc = Eigen::MatrixXd::Identity(2,2);
    R =  Eigen::MatrixXd::Identity(2,2);
    M = Eigen::MatrixXd::Zero(2,2); //0.25 * Eigen::MatrixXd::Identity(2,2);
    kmath::LinearSystem sys(Ac, Bc, Eigen::MatrixXd());
    std::cout << "Controllable: " << sys.is_controllable() << "\n";

    std::cout <<"-----------------------LQR test------------------------\n";
    start = kite_utils::get_time();
    Eigen::MatrixXd K = kmath::oc::lqr(sys, Qc, R, M);
    finish = kite_utils::get_time();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    latency = static_cast<double>(duration.count());
    std::cout << "Feedback matrix: \n" << K << "\n";
    std::cout << "Time elapsed : " << latency * 1e-6 << " [seconds]" << "\n";

    BOOST_CHECK(true);
}

SX ode(const SX &x, const SX &u, const SX &p)
{
    SX f = SX::zeros(2,1);
    f(0) = x(1) + u(0);
    f(1) = u(1);

    return f;
}

BOOST_AUTO_TEST_CASE( pseudo_test )
{
    const int poly_order   = 2;
    const int num_segments = 3;
    const int dimx         = 2;
    const int dimu         = 2;
    const int dimp         = 0;
    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> cheb;
    std::cout << "----------------------------------------------- \n";
    std::cout << "Collocation Matrix \n" << cheb.CPoints() << " \n";
    std::cout << "Differentiation Matrix \n" << cheb.D() << "\n";
    std::cout << "Quadrature weights \n" << cheb.QWeights() << "\n";
    std::cout << "Composite Diff matrix \n" << cheb.CompD() << "\n";
    std::cout << "----------------------------------------------- \n";

    SX x = SX::sym("x", 2);
    SX u = SX::sym("u", 2);
    SX f = SX::vertcat({x(1) + u(0), u(1)});
    std::cout << "F : " << f << "\n";
    Function dynamics = Function("rhs", {x, u}, {f});
    std::cout << "X var : " << cheb.VarX() << "\n";
    std::cout << "U var : " << cheb.VarU() << "\n";


    SX g = cheb.CollocateDynamics(dynamics, 0, 1);
    auto G  = cheb.CollocateDynamics2(ode, 0, 1);
    SX g_fun = G(cheb.VarX(), cheb.VarU(), SX({0}));


    std::cout << "Collocated Dynamics \n" << g << "\n";
    std::cout << "Collocated Dunamcis2 \n" << g_fun << "\n";

    /** Performance index dicretization */
    SX sym_M = SX::norm_2(x);
    SX sym_L = SX::norm_2(x) + SX::norm_2(u);

    Function M = Function(); //Function("Mayer", {x}, {sym_M});
    Function L = Function("Lagrange", {x, u}, {sym_L});
    SX J = cheb.CollocateCost(M, L, 0, 1);
    std::cout << "Collocated Cost :  \n" << J << "\n";
    std::cout << "Cost size: " << J.size() << "\n";

    std::cout << "Chebyshev Expansion  : \n";
    std::vector<double> coeff = {1, 2, 3, 4};
    double value = -1;
    double expansion = kmath::chebyshev_expansion<double>(coeff, value);
    double expansion2 = kmath::chebyshev_expansion2<double>(coeff, value);
    std::cout << expansion << "\n";
    std::cout << expansion2 << "\n";

    DM cp, D;
    kmath::cheb(cp, D, poly_order, std::make_pair<double>(-1,1));
    std::cout << "ref collocation points: \n" << cp << "\n";
    std::cout << "ref diff matrix: \n"        << D << "\n";

    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE( collocation_test )
{
    /** load data from the log file */
    /** std::ifstream log_file;
    int nrows = 100;
    int ncols = 14;
    std::vector<std::vector<double>> data_matrix(nrows, std::vector<double>(ncols));
    log_file.open("kite_state_cut.log", std::ios::in);
    if(!log_file)
        std::cerr << "Unable to open lof file: " << "kite_state_cut.log" << "\n";

    for (int i = 0; i < nrows; ++i)
        for(int j = 0; j < ncols; ++j)
            log_file >> data_matrix[i][j];

    DM kite_states(data_matrix);

    /** create a controller and solve an OCP */
    std::string kite_config_file = "umx_radian.yaml";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_config_file);
    AlgorithmProperties algo_props;
    algo_props.Integrator = RK4;

    /** set up the path */
    SX x = SX::sym("x");
    double radius   = 2.31;
    double altitude = 0.00;
    SX Path = SX::vertcat(SXVector{radius * cos(x), radius * sin(x), altitude});
    Function path_fun = Function("path", {x}, {Path});
    std::shared_ptr<KiteDynamics> kite = std::make_shared<KiteDynamics>(kite_props, algo_props);

    /** Kite Dynamics debugging */
    Function DYNAMO = kite->getNumericDynamics();
    DM x0 = DM({1.5, 0, 0, 0, 0, 0, 0, 1.0, 0, 1, 0, 0.0, 0.0});
    DM u0 = DM({0.1, 0, 0});

    DMVector eval = DYNAMO(DMVector{x0,u0});
    std::cout << "Kite Dynamics: " << eval[0] << "\n";

    KiteNMPF controller = KiteNMPF(kite, path_fun);

    /** set control constraints */
    double angle_sat = kmath::deg2rad(8.0);
    DM lbu = DM::vertcat({0, -angle_sat, -angle_sat, -10});
    DM ubu = DM::vertcat({0.3, angle_sat, angle_sat, 10});
    controller.setLBU(lbu);
    controller.setUBU(ubu);

    /** set variable constraints */
    DM lbx = DM::vertcat({0.5, -0.5, -DM::inf(1), -2 * M_PI, -2 * M_PI, -2 * M_PI, -DM::inf(1), -DM::inf(1), -DM::inf(1),
                          -1, -1, -1, -1, -DM::inf(1), -DM::inf(1)});

    DM ubx = DM::vertcat({12, 5, DM::inf(1), 2 * M_PI, 2 * M_PI, 2 * M_PI, DM::inf(1), DM::inf(1), DM::inf(1),
                          1, 1, 1, 1, DM::inf(1), DM::inf(1)});
    controller.setLBX(lbx);
    controller.setUBX(ubx);

    /** set reference velocity */
    DM vel_ref = 0.05;
    controller.setReferenceVelocity(vel_ref);

    /** set scaling of system ODEs */
    double vx_max = ubx.nonzeros()[0];
    double vy_max = ubx.nonzeros()[1];
    double wx_max = ubx.nonzeros()[3];
    double wy_max = ubx.nonzeros()[4];
    double wz_max = ubx.nonzeros()[5];
    double rsphere = 5;

    //DM Sx = DM::diag(DM::vertcat(DMVector{1/vx_max, 1/vy_max, 1, 1/wx_max, 1/wy_max, 1/wz_max,
    //                            1/rsphere, 1/rsphere, 1/rsphere, 1, 1, 1, 1, 1/(2*M_PI), 1/vx_max}));
    DM Sx = DM::eye(15);
    /** control */
    DM Su = DM::eye(4);
    controller.setStateScaling(Sx);
    controller.setControlScaling(Su);

    /** create path following NLP */
    controller.createNLP();

    /** solve for initial conditions */
    DM X0 = DM::vertcat({1.5, 0, 0, 0, 0, 0,
                         0, -1.0, 0, 1, 0, 0.0, 0.0, M_PI_2, 0});

    std::chrono::time_point<std::chrono::system_clock> start = kite_utils::get_time();

    controller.computeControl(X0);

    std::chrono::time_point<std::chrono::system_clock> stop = kite_utils::get_time();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    auto latency = static_cast<double>(duration.count());

    DM opt_x = controller.getOptimalTrajetory();
    std::cout << "Optimal Solution: \n" << opt_x << "\n";

    DM opt_u = controller.getOptimalControl();
    std::cout << "Optimal Control : \n" << opt_u << "\n";
    std::cout << "Actual Control : \n" << opt_u(Slice(0, 4), opt_u.size2() - 1) << "\n";

    DM U0 = DM::zeros(4,1);
    for(int i = 0; i < opt_u.size2()-1; ++i)
        U0 += pow(-1, i) * opt_u(Slice(0,4), i);
    std::cout << "Extrapolated control : \n" << U0 << "\n";

    std::cout << "COLLOCATION_TEST COMPUTATION TIME: " << latency * 1e-6 << " [seconds]" << "\n";

    /** initialization */
    DM theta0 = controller.findClosestPointOnPath(X0(Slice(6,9)));
    std::cout << "INIT: " << theta0 << "\n";

    /** process real trajectory */
    bool success = true;
    //for(int i = 0; i < kite_states.size1(); ++i)
    /** for(int i = 0; i < 1; ++i)
    {
        start = kite_utils::get_time();
        X0 = kite_states(i, Slice(1, kite_states.size2()));
        X0 = DM::horzcat({X0, 0, 0}).T();
        controller.computeControl(X0);
        stop = kite_utils::get_time();

        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        latency = static_cast<double>(duration.count());
        std::cout << "COLLOCATION_TEST COMPUTATION TIME: " << latency * 1e-6 << " [seconds]" << "\n";
    } */

    BOOST_CHECK(success);
}


BOOST_AUTO_TEST_CASE( simple_ocp_test )
{
    const int poly_order   = 3;
    const int num_segments = 1;
    const int dimx         = 1;
    const int dimu         = 1;
    const int dimp         = 0;
    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> cheb;

    SX x = SX::sym("x");
    SX u = SX::sym("u");
    SX f = -x + u;
    Function dynamics = Function("rhs", {x, u}, {f});
    SX g = cheb.CollocateDynamics(dynamics, 0, 1);

    /** Performance index dicretization */
    SX sym_M = x;
    SX sym_L = pow(u, 4);

    Function M = Function("Mayer", {x}, {sym_M});
    Function L = Function("Lagrange", {x, u}, {sym_L});
    SX J = cheb.CollocateCost(M, L, 0, 1);

    SX varx = cheb.VarX();
    SX varu = cheb.VarU();

    SX opt_var = SX::vertcat(SXVector{varx, varu});

    SX lbg = SX::zeros(g.size());
    SX ubg = SX::zeros(g.size());

    /** set inequality (box) constraints */
    /** state */
    SX lbx = SX::repmat(-SX::inf(), poly_order + 1, 1);
    SX ubx = SX::repmat(SX::inf(), poly_order + 1, 1);

    lbx = SX::vertcat( {lbx, SX::repmat(-SX::inf(), poly_order + 1, 1)} );
    ubx = SX::vertcat( {ubx, SX::repmat(SX::inf(), poly_order + 1, 1)} );

    /** formulate NLP */
    SXDict NLP;
    NLP["x"] = opt_var;
    NLP["f"] = J;
    NLP["g"] = g;

    Dict OPTS;
    OPTS["ipopt.linear_solver"]  = "ma97";
    OPTS["ipopt.print_level"]    = 0;
    OPTS["ipopt.tol"]            = 1e-5;
    OPTS["ipopt.acceptable_tol"] = 1e-4;
    //OPTS["ipopt.max_iter"]       = 15;

    Function NLP_Solver;
    NLP_Solver = nlpsol("solver", "ipopt", NLP, OPTS);

    /** set default args */
    DMDict ARG;
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    /** solve */
    DMDict res = NLP_Solver(ARG);
    std::cout << NLP_Solver.stats() << "\n";

    BOOST_CHECK(true);
}

SX add(const SX &a, const SX &b, const SX &c)
{
    return a + b;
}

BOOST_AUTO_TEST_CASE( simple_generics_test )
{
    const int poly_order   = 3;
    const int num_segments = 1;
    const int dimx         = 1;
    const int dimu         = 1;
    const int dimp         = 0;
    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> cheb;

    SX x = SX::sym("x");
    SX u = SX::sym("u");
    SX p = SX::sym("p");

    auto f3 = cheb.CollocateDynamics2(add, 0, 1);
    std::cout << "F2 : " << cheb._ode(x, u, p) << "\n";
    //std::cout << "F3 : " << f3(x, u, p) << "\n";

    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_CASE( full_generics_test )
{
    /** define kite dynamics */
    std::string kite_params_file = "umx_radian.yaml";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_params_file);
    AlgorithmProperties algo_props;
    algo_props.Integrator = RK4;
    algo_props.sampling_time = 0.02;
    std::shared_ptr<KiteDynamics> kite = std::make_shared<KiteDynamics>(kite_props, algo_props);

    /** define path function */
    SX y = SX::sym("y");
    double radius   = 3.0;
    double altitude = 0.00;
    SX Path = SX::vertcat(SXVector{radius * cos(y), radius * sin(y), altitude});
    /** rotate path */
    SX q_rot = SX::vertcat({cos(M_PI / 24), 0, sin(M_PI / 24), 0});
    SX q_rot_inv = kmath::quat_inverse(q_rot);
    SX qP_tmp = kmath::quat_multiply(q_rot_inv, SX::vertcat({0, Path}));
    SX qP_q = kmath::quat_multiply(qP_tmp, q_rot);
    Path = qP_q(Slice(1,4), 0);
    Function PathFunc = Function("path", {y}, {Path});

    /**-------------------------------------------------------------------*/

    /** get dynamics function and state Jacobian */
    SX dynamics = kite->getSymbolicDynamics();
    SX X = kite->getSymbolicState();
    SX U = kite->getSymbolicControl();

    /** state and control dimensionality */
    int n = 15;
    int m = 4;

    /** Order of polynomial interpolation */
    int N = 10;

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

    /** evaluate augmented dynamics */
    Function aug_dynamo = Function("AUG_DYNAMO", {aug_state, aug_control}, {aug_dynamics});
    Function DynamicsFunc = aug_dynamo;

    DM Q  = 1e-2 * SX::diag(SX({1e4, 1e4, 5e3}));
    DM R  = 1e-4 * SX::diag(SX({1, 1, 1, 1}));
    DM W  = 1e-3;

    DM LBX = -DM::inf(15);
    DM UBX = DM::inf(15);

    DM LBU = -DM::inf(4);
    DM UBU = DM::inf(4);

    SX reference_velocity = 0.05;

    /** ----------------------------------------------------------------------------------*/
    const int num_segments = 1;
    const int poly_order   = 10;
    const int dimx         = 15;
    const int dimu         = 4;
    const int dimp         = 0;

    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> spectral;
    SX diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, 1);

    /** define an integral cost */
    SX x = SX::sym("x", dimx);
    SX u = SX::sym("u", dimu);

    SXVector tmp = PathFunc(SXVector{x(13)});
    SX sym_path  = tmp[0];
    SX residual  = sym_path - x(Slice(6,9));
    SX lagrange  = SX::sum1( SX::mtimes(Q, pow(residual, 2)) ) + SX::sum1( SX::mtimes(W, pow(reference_velocity - x(14), 2)) );
    Function LagrangeTerm = Function("Lagrange", {x, u}, {lagrange});

    SX mayer     =  SX::sum1( SX::mtimes(2 * Q, pow(residual, 2)) );
    Function MayerTerm    = Function("Mayer",{x}, {mayer});
    SX performance_idx = spectral.CollocateCost(MayerTerm, LagrangeTerm, 0, 1);

    SX varx = spectral.VarX();
    SX varu = spectral.VarU();

    SX opt_var = SX::vertcat(SXVector{varx, varu});

    SX lbg = SX::zeros(diff_constr.size());
    SX ubg = SX::zeros(diff_constr.size());

    /** set inequality (box) constraints */
    /** state */
    SX lbx = SX::repmat(LBX, poly_order + 1, 1);
    SX ubx = SX::repmat(UBX, poly_order + 1, 1);

    /** control */
    lbx = SX::vertcat( {lbx, SX::repmat(LBU, poly_order + 1, 1)} );
    ubx = SX::vertcat( {ubx, SX::repmat(UBU, poly_order + 1, 1)} );

    SX diff_constr_jacobian = SX::jacobian(diff_constr, opt_var);
    /** Augmented Jacobian */
    Function AugJacobian = Function("aug_jacobian",{opt_var}, {diff_constr_jacobian});

    /** formulate NLP */
    SXDict NLP;
    Dict OPTS;
    DMDict ARG;
    NLP["x"] = opt_var;
    NLP["f"] = performance_idx;
    NLP["g"] = diff_constr;

    OPTS["ipopt.linear_solver"]  = "ma97";
    OPTS["ipopt.print_level"]    = 0;
    OPTS["ipopt.tol"]            = 1e-5;
    OPTS["ipopt.acceptable_tol"] = 1e-4;
    OPTS["ipopt.max_iter"]       = 20;
    Function NLP_Solver = nlpsol("solver", "ipopt", NLP, OPTS);

    /** set default args */
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    DM feasible_state =  (UBX + LBX) / 2;
    DM feasible_control = (UBU + LBU) / 2;

    //ARG["x0"] = DM::vertcat(DMVector{DM::repmat(feasible_state, poly_order + 1, 1),
    //                                 DM::repmat(feasible_control, poly_order + 1, 1)});
    ARG["x0"] = DM::vertcat({0.322159, -1.7086, 0.0187737, 0.571611, -0.463085, -3.98942, -1.00108, 2.67555, -1.67471, 0.685529, 0.179051, 0.601294,
            -0.234078, -8.93944, -6.15368, 0.320558, -1.81794, -0.048227, 0.697989, -0.487516, -3.89813, -0.978354, 2.70897, -1.6673, 0.693579,
            0.203702, 0.597764, -0.194959, -8.78899, -6.13793, 0.254425, -1.85413, 0.141461, 0.852417, -0.100141, -3.66057, -0.922693, 2.80918,
            -1.6289, 0.712096, 0.26544, 0.576425, -0.0820664, -8.3567, -6.02738, 7.80537e-09, -2.46967, 2.21734e-09, 1.44569, -0.464535,
            -2.70741, -0.862032, 2.98398, -1.51414, 0.725987, 0.318053, 0.529904, 0.0847142, -7.7014, -5.81885, 0.143081, -2.19291, 0.185359,
            1.49807, -1.06445, -1.31733, -0.826037, 3.2155, -1.30218, 0.705787, 0.300448, 0.509612, 0.259055, -6.9092, -5.54832, 0.455407,
            -1.48802, 0.925243, 0.939966, -0.431242, -0.183387, -0.971463, 3.33088, -1.10238, 0.683697, 0.242655, 0.507753, 0.363209, -6.07502,
            -5.24718, 1.30346, -1.11843, 0.985835, 0.52808, 0.239268, 0.253541, -1.13998, 3.37537, -0.955375, 0.696112, 0.198938, 0.49205,
            0.386145, -5.28823, -4.94104, 3.23749, -0.562678, 1.1487, 0.185406, 0.798954, 0.427541, -1.37814, 3.14841, -0.761349, 0.725335,
            0.182191, 0.46157, 0.378231, -4.61785, -4.66725, 4.98206, 0.581094, 1.0498, -0.361074, 1.36802, -0.142782, -1.43793, 2.9391,
            -0.593935, 0.758427, 0.200775, 0.420269, 0.347893, -4.11561, -4.44789, 5.96975, 0.969053, 0.86408, -0.503073, 1.45039, -0.880998,
            -1.73889, 2.39822, -0.451459, 0.761644, 0.241508, 0.383887, 0.382548, -3.79856, -4.30918, 6.17882, 1.08634, 0.852976, -0.562302,
            1.54206, -1.17198, -0.95457, 2.91599, -0.201047, 0.867553, 0.374954, 0.230215, 0.231888, -3.79856, -4.30918, 0.121739, 0.136834,
            0.136834, -0.165329, 0.101, 0.116373, -0.136834, -1.03451, 0.199, -0.0223408, 0.0984078, -1.81909, 0.129195, -0.136834, 0.136834,
            -1.91874, 0.101, -0.135728, 0.0148018, -1.94518, 0.101, -0.106827, 0.0835699, -1.96381, 0.122227, -0.136834, -0.115507, -1.9685,
            0.101, -0.136834, -0.136834, -1.96848, 0.199, -0.135985, -0.128318, -1.9641, 0.101, -0.136834, 0.00392875, -1.93056, 0.199,
            -0.136834, 0.0482633, -1.83182});

    std::cout << "Try to detect corrupted jacobain: \n";
    DMVector res  = AugJacobian(ARG["x0"]);
    std::cout << res[0] << "\n";

    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_SUITE_END()

