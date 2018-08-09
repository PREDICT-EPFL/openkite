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
    std::string boat_config_file = "config_id.yaml";
    BoatProperties boat_props = BoatProperties::Load(boat_config_file);

    bifoiler::BoatDynamics boat_int(boat_props);    //integration model
    Function ode = boat_int.getNumericDynamics();

    /** compare three ode solvers */
    Dict opts;
    opts["tf"]         = 3.0;
    opts["poly_order"] = 100;
    opts["tol"]        = 1e-4;
    opts["method"]     = IntType::RK4;
    ODESolver rk4_solver(ode, opts);

    opts["method"] = IntType::CVODES;
    ODESolver cvodes_solver(ode, opts);

    opts["method"] = IntType::CHEBYCHEV;
    ODESolver chebychev_solver(ode, opts);

    double tf = 7.0;
    casadi::DMDict props;
    props["scale"] = 0;
    props["P"] = casadi::DM::diag(casadi::DM({0.1, 1/3.0, 1/3.0, 1/2.0, 1/5.0, 1/2.0, 1/3.0, 1/3.0, 1/3.0, 1.0, 1.0, 1.0, 1.0}));
    props["R"] = casadi::DM::diag(casadi::DM({1/0.3, 1/0.3, 1/0.3, 1}));
    PSODESolver<20,20,13,4>ps_solver(ode, tf, props);

    std::cout << "All solvers have been created \n";

    /** solve a problem */
    DM rk4_sol, cheb_sol, cv_sol, ps_sol;

    DM init_state = DM::vertcat({5.4663474e+00,   8.6672711e-01,   6.5021107e-02,  -1.0385855e-01,  -6.0455710e-02 ,  2.0441055e-01,
                                 -2.3351574e+00,   4.3712858e+01,   9.3579102e-01,  -3.0662162e-02,   3.0009789e-02,  -7.5182736e-01,  -6.5796280e-01});
    DM control = DM::vertcat({-4.8015300e-02,   1.2532470e-01,   3.0563000e-02,   6.8800000e-01});

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
    std::cout << "CVODES solve time: " << std::setprecision(6)
              << static_cast<double>(cv_duration.count()) * 1e-6 << " [seconds]" << "\n";
    std::cout << "CHEB solve time: " << std::setprecision(6)
              << static_cast<double>(cheb_duration.count()) * 1e-6 << " [seconds]" << "\n";
    std::cout << "PS solve time: " << std::setprecision(6)
              << static_cast<double>(ps_duration.count()) * 1e-6 << " [seconds]" << "\n";

    std::cout << "RK4: " <<  rk4_sol << "\n";
    //std::cout << "CVODES: " <<  cv_sol << "\n";
    std::cout << "CHEB: " <<  cheb_sol << "\n";
    std::cout << "PS:" << ps_sol(Slice(0, 13)) << "\n";

    std::cout << "SIZE: " << ps_sol.size1() << "\n";

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
