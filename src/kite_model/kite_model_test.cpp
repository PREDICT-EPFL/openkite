#define BOOST_TEST_TOOLS_UNDER_DEBUGGER
#define BOOST_TEST_MODULE kite_model_test
#include <boost/test/included/unit_test.hpp>

#include "integrator.h"
#include "kite.h"

using namespace casadi;

BOOST_AUTO_TEST_SUITE( kite_model_suite_test )

BOOST_AUTO_TEST_CASE( ode_solver_test )
{
    std::string kite_config_file = "umx_radian.yaml";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_config_file);
    AlgorithmProperties algo_props;
    algo_props.Integrator = RK4;

    /** SCALING */
    //kite_props.Tether.length = kite_props.Tether.length / 3;

    KiteDynamics kite(kite_props, algo_props);
    Function ode = kite.getNumericDynamics();

    /** compare three ode solvers */
    Dict opts;
    opts["tf"]         = 5.0;
    opts["poly_order"] = 41;
    opts["tol"]        = 1e-4;
    opts["method"]     = IntType::RK4;
    ODESolver rk4_solver(ode, opts);

    opts["method"] = IntType::CVODES;
    ODESolver cvodes_solver(ode, opts);

    opts["method"] = IntType::CHEBYCHEV;
    ODESolver chebychev_solver(ode, opts);


    //kite_props.Tether.length = kite_props.Tether.length / 3;
    //std::cout << "New tether length : " << kite_props.Tether.length << "\n";

    KiteDynamics kite2(kite_props, algo_props);
    Function ode2 = kite.getNumericDynamics();
    double tf = 7.0;
    casadi::DMDict props;
    props["scale"] = 0;
    props["P"] = casadi::DM::diag(casadi::DM({0.1, 1/3.0, 1/3.0, 1/2.0, 1/5.0, 1/2.0, 1/3.0, 1/3.0, 1/3.0, 1.0, 1.0, 1.0, 1.0}));
    props["R"] = casadi::DM::diag(casadi::DM({1/0.15, 1/0.2618, 1/0.2618}));
    PSODESolver<10,40,13,3>ps_solver(ode, tf, props);

    /** solve a problem */
    DM rk4_sol, cheb_sol, cv_sol, ps_sol;
    //DM init_state = DM::vertcat({0.318732, 0.182552, 0.254833, 1.85435, -0.142882, -0.168359,
    //                          -0.229383, -0.0500282, -0.746832, 0.189409, -0.836349, -0.48178, 0.180367});
    //DM control = DM::vertcat({0.3, kmath::deg2rad(5), -kmath::deg2rad(2)});

    DM init_state = DM::vertcat({6.1977743e+00,  -2.8407148e-02,   9.1815942e-01,   2.9763089e-01,  -2.2052198e+00,  -1.4827499e-01,
                                 -4.1624807e-01, -2.2601052e+00,   1.2903439e+00,   3.5646195e-02,  -6.9986094e-02,   8.2660637e-01,   5.5727089e-01});
    DM control = DM::vertcat({0.1, 0.0, 0.0});

    std::chrono::time_point<std::chrono::system_clock> start = kite_utils::get_time();

    rk4_sol  = rk4_solver.solve(init_state, control, tf);
    std::chrono::time_point<std::chrono::system_clock> rk4_stop = kite_utils::get_time();

    //cv_sol   = cvodes_solver.solve(init_state, control, tf);
    std::chrono::time_point<std::chrono::system_clock> cv_stop = kite_utils::get_time();

    cheb_sol = chebychev_solver.solve(init_state, control, tf);
    std::chrono::time_point<std::chrono::system_clock> cheb_stop = kite_utils::get_time();

    bool FULL = true;
    ps_sol = ps_solver.solve(init_state, control, FULL);
    std::chrono::time_point<std::chrono::system_clock> ps_stop = kite_utils::get_time();

    auto rk4_duration = std::chrono::duration_cast<std::chrono::microseconds>(rk4_stop - start);
    auto cv_duration = std::chrono::duration_cast<std::chrono::microseconds>(cv_stop - rk4_stop);
    auto cheb_duration = std::chrono::duration_cast<std::chrono::microseconds>(cheb_stop - cv_stop);
    auto ps_duration = std::chrono::duration_cast<std::chrono::microseconds>(ps_stop - cheb_stop);

    std::cout << "RK4 solve time: " << std::setprecision(6)
              << static_cast<double>(rk4_duration.count()) * 1e-6 << " [seconds]" << "\n";
    //std::cout << "CVODES solve time: " << std::setprecision(6)
    //          << static_cast<double>(cv_duration.count()) * 1e-6 << " [seconds]" << "\n";
    std::cout << "CHEB solve time: " << std::setprecision(6)
              << static_cast<double>(cheb_duration.count()) * 1e-6 << " [seconds]" << "\n";
    std::cout << "PS solve time: " << std::setprecision(6)
              << static_cast<double>(ps_duration.count()) * 1e-6 << " [seconds]" << "\n";

    std::cout << "RK4: " <<  rk4_sol << "\n";
    std::cout << "CVODES: " <<  cv_sol << "\n";
    std::cout << "CHEB: " <<  cheb_sol << "\n";
    std::cout << "PS:" << ps_sol(Slice(0, 13)) << "\n";

    std::ofstream trajectory_file("integrated_trajectory.txt", std::ios::out);

    if(!trajectory_file.fail())
    {
        for (int i = 0; i < ps_sol.size1(); i = i + 13)
        {
            std::vector<double> tmp = ps_sol(Slice(i, i + 13),0).nonzeros();
            for (uint j = 0; j < tmp.size(); j++)
            {
                trajectory_file << tmp[j] << " ";
            }
            trajectory_file << "\n";
        }
    }
    trajectory_file.close();

    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
