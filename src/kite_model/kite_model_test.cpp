#define BOOST_TEST_TOOLS_UNDER_DEBUGGER
#define BOOST_TEST_MODULE kite_model_test
#include <boost/test/included/unit_test.hpp>

#include "integrator.h"
#include "boat_model.h"
#include "kite.h"

using namespace casadi;

BOOST_AUTO_TEST_SUITE( kite_model_suite_test )

BOOST_AUTO_TEST_CASE( ode_solver_test )
{   
    std::string boat_config_file = "config.yaml";
    BoatProperties boat_props = BoatProperties::Load(boat_config_file);

    bifoiler::BoatDynamics boat(boat_props, true);
    bifoiler::BoatDynamics boat_int(boat_props);    //integration model
    Function ode = boat_int.getNumericDynamics();

    /** compare three ode solvers */
    Dict opts;
    opts["tf"]         = 1.0;
    opts["poly_order"] = 21;
    opts["tol"]        = 1e-4;
    opts["method"]     = IntType::RK4;
    ODESolver rk4_solver(ode, opts);

    opts["method"] = IntType::CVODES;
    ODESolver cvodes_solver(ode, opts);

    opts["method"] = IntType::CHEBYCHEV;
    ODESolver chebychev_solver(ode, opts);

    double tf = 1.0;
    casadi::DMDict props;
    props["scale"] = 0;
    props["P"] = casadi::DM::diag(casadi::DM({0.1, 1/3.0, 1/3.0, 1/2.0, 1/5.0, 1/2.0, 1/3.0, 1/3.0, 1/3.0, 1.0, 1.0, 1.0, 1.0}));
    props["R"] = casadi::DM::diag(casadi::DM({1/0.15, 1/0.2618, 1/0.2618}));
    PSODESolver<10,4,13,3>ps_solver(ode, tf, props);

    /** solve a problem */
    DM rk4_sol, cheb_sol, cv_sol, ps_sol;

    DM init_state = DM::vertcat({7.2547334e+00, -5.9953832e-01, 1.5621571e+00, 6.4458230e-01, -1.9436366e+00,
                                 -1.9995872e+00, 1.6915701e+00, -2.0761443e+00, -3.6169709e-01, 4.0336430e-01, 1.5472226e-01, -3.7687924e-01, -8.1934426e-01});
    DM control = DM::vertcat({0.15, 0.0, 0.0});

    std::chrono::time_point<std::chrono::system_clock> start = kite_utils::get_time();

    rk4_sol  = rk4_solver.solve(init_state, control, tf);
    std::chrono::time_point<std::chrono::system_clock> rk4_stop = kite_utils::get_time();

    cv_sol   = cvodes_solver.solve(init_state, control, tf);
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
    std::cout << "CVODES solve time: " << std::setprecision(6)
              << static_cast<double>(cv_duration.count()) * 1e-6 << " [seconds]" << "\n";
    std::cout << "CHEB solve time: " << std::setprecision(6)
              << static_cast<double>(cheb_duration.count()) * 1e-6 << " [seconds]" << "\n";
    std::cout << "PS solve time: " << std::setprecision(6)
              << static_cast<double>(ps_duration.count()) * 1e-6 << " [seconds]" << "\n";

    std::cout << "RK4: " <<  rk4_sol << "\n";
    std::cout << "CVODES: " <<  cv_sol << "\n";
    std::cout << "CHEB: " <<  cheb_sol << "\n";
    std::cout << "PS:" << ps_sol(Slice(0, 13)) << "\n";

    std::ofstream trajectory_file("estimated_trajectory.txt", std::ios::out);

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
