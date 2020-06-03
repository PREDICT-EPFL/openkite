#define BOOST_TEST_TOOLS_UNDER_DEBUGGER
#define BOOST_TEST_MODULE kite_model_test

#include <boost/test/included/unit_test.hpp>

#include "integrator.h"
#include "kite.h"

using namespace casadi;

BOOST_AUTO_TEST_SUITE(kite_model_suite_test)

    BOOST_AUTO_TEST_CASE(ode_solver_test) {

        std::string flightDataDir = "/home/johannes/identification/processed_flight_data/";

        struct FlightMetaData {
            /// 1. ///
            int session = 2;
            int number = 1;
            int seq = 1;

            /// 2. ///
            std::string maneuver = "identRoll";
            //std::string maneuver = "identPitch";

            /// 3. ///
            // (always equal!) std::string resampleMethod = "cheb";
            std::string resampleMethod = "equal";
        } flight;

        /// 4. ///
        const int DATA_POINTS = 90;

        std::string flightDataPath = flightDataDir
                                     + "session_" + std::to_string(flight.session)
                                     + "/flight_" + std::to_string(flight.number);
        if (!flight.maneuver.empty()) flightDataPath.append("/" + flight.maneuver);
        flightDataPath.append("/seq_" + std::to_string(flight.seq) + "/");

        std::cout << flightDataPath << "\n";

        /** define kite dynamics */
        /// 5. ///
        std::string kite_config_file;
        kite_config_file = "/home/johannes/identification/eg4.yaml";
        KiteProperties kite_props = kite_utils::LoadProperties(kite_config_file);

        /** Load wind data **/
        std::ifstream id_wind_file(flightDataPath + "wind.txt", std::ios::in);
        if (!id_wind_file.fail()) {
            int windFrom_deg;
            double windSpeed;
            id_wind_file >> windFrom_deg;
            id_wind_file >> windSpeed;

            kite_props.atmosphere.WindFrom = windFrom_deg * M_PI/180.0;
            kite_props.atmosphere.WindSpeed = windSpeed;
            kite_props.atmosphere.airDensity= 1.1589; // Standard atmosphere at 468 meters
        } else {
            std::cout << "Could not open : id wind data file \n";
            id_wind_file.clear();
        }

        /** Load validation data */
        std::ifstream id_data_file(flightDataPath + flight.resampleMethod + "_states.txt", std::ios::in);
        std::ifstream id_control_file(flightDataPath + flight.resampleMethod + "_controls.txt", std::ios::in);
        const int state_size = 13;  // v(3) w(3) r(3) q(4)
        const int control_size = 4; // T elev rud ail

        DM id_data = DM::zeros(1 + state_size, DATA_POINTS);        // Time + States
        DM id_control = DM::zeros(1 + control_size, DATA_POINTS);   // Time + Controls

        /** load state trajectory */
        if (!id_data_file.fail()) {
            for (uint i = 0; i < DATA_POINTS; ++i) {
                for (uint j = 0; j < (1 + state_size); ++j) {       // Time + States
                    double entry;
                    id_data_file >> entry;
                    id_data(j, i) = entry;
                }
            }
        } else {
            std::cout << "Could not open : id state data file \n";
            id_data_file.clear();
        }

        /** load control data */
        if (!id_control_file.fail()) {
            for (uint i = 0; i < DATA_POINTS; ++i) {
                for (uint j = 0; j < (1 + control_size); ++j) {     // Time + Controls
                    double entry;
                    id_control_file >> entry;
                    /** put in reverse order to comply with Chebyshev method */
                    id_control(j, DATA_POINTS - 1 - i) = entry;
                }
            }
        } else {
            std::cout << "Could not open : id control data file \n";
            id_control_file.clear();
        }

        AlgorithmProperties algo_props;
        algo_props.Integrator = RK4;

        // End time - start time of data
        const double timeStep_controls = static_cast<double>(id_control(0, id_control.size2() - 2) -
                                                             id_control(0, id_control.size2() - 1));
        std::cout << "timeStep_controls: " << timeStep_controls << "\n";
        //std::cout << id_control.size() << "\n";

        /** SCALING */
        //kite_props.Tether.length = kite_props.Tether.length / 3;

        KiteDynamics kite(kite_props, algo_props);
        Function ode = kite.getNumericDynamics();

        /** compare three ode solvers */
        Dict opts;
        opts["tf"] = timeStep_controls; //1.0;
        opts["poly_order"] = 41;
        opts["tol"] = 1e-4;
        opts["method"] = IntType::RK4;
        ODESolver rk4_solver(ode, opts);

        opts["method"] = IntType::CVODES;
        ODESolver cvodes_solver(ode, opts);

        opts["method"] = IntType::CHEBYCHEV;
        ODESolver chebychev_solver(ode, opts);


        //kite_props.Tether.length = kite_props.Tether.length / 3;
        //std::cout << "New tether length : " << kite_props.Tether.length << "\n";

        KiteDynamics kite2(kite_props, algo_props);
        Function ode2 = kite.getNumericDynamics();
        double tf = timeStep_controls; //1.0;
        casadi::DMDict props;
        props["scale"] = 0;
        props["P"] = casadi::DM::diag(casadi::DM(
                {0.1, 1 / 3.0, 1 / 3.0,
                 1 / 2.0, 1 / 5.0, 1 / 2.0,
                 1 / 3.0, 1 / 3.0, 1 / 3.0,
                 1.0, 1.0, 1.0, 1.0}));
        props["R"] = casadi::DM::diag(casadi::DM({1 / 0.15, 1 / 0.2618, 1 / 0.2618, 1 / 0.2618}));
        PSODESolver<3, 42, 13, 4> ps_solver(ode, tf, props);

        //            DM init_state = DM::vertcat({15, -0, 0, // vx vy vz
        //                                         0, 0, 0, // wx wy wz
        //                                         0, 0, -2, // x y z
        //                                         3.5646195e-02, -6.9986094e-02, 8.2660637e-01, 5.5727089e-01 // qw qx qy qz
        //                                        });
        DM init_state = id_data(Slice(1, id_data.size1()), 0);
        std::cout << init_state << "\n";

        DM sim_states = DM::zeros(1 + state_size, DATA_POINTS);         // Time + States x Data points (discrete times)
        sim_states(0, Slice(0, sim_states.size2())) = id_data(0, Slice(0, sim_states.size2())); // Copy times
        sim_states(Slice(1, 1 + state_size), 0) = id_data(Slice(1, id_data.size1()),
                                                          0); // Copy initial state (w/o time)

        for (int iDatapoint = 0; iDatapoint < (id_data.size2() - 1); ++iDatapoint) {
            /** solve a problem */
            DM rk4_sol, cheb_sol, cv_sol, ps_sol;

            std::cout << sim_states << "\n";
            DM state = sim_states(Slice(1, sim_states.size1()), iDatapoint);

            //            DM control = DM::vertcat({0.1, 0.0, 0.0, 0.0});
            //DM control = id_control(Slice(1, id_control.size1()), iDatapoint);
            DM control = id_control(Slice(1, id_control.size1()), DATA_POINTS - 1 - iDatapoint);

            std::cout << "state: " << state << "\n";
            std::cout << "control: " << control << "\n";

            std::chrono::time_point<std::chrono::system_clock> start = kite_utils::get_time();
            //TODO
            rk4_sol = rk4_solver.solve(state, control, tf);
            std::chrono::time_point<std::chrono::system_clock> rk4_stop = kite_utils::get_time();

            //TODO
            cv_sol = cvodes_solver.solve(init_state, control, tf);
            std::chrono::time_point<std::chrono::system_clock> cv_stop = kite_utils::get_time();

            //TODO
            //cheb_sol = chebychev_solver.solve(state, control, tf);
            std::chrono::time_point<std::chrono::system_clock> cheb_stop = kite_utils::get_time();

            bool FULL = true;
            ps_sol = ps_solver.solve(state, control, FULL);
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

            std::cout << "RK4: " << rk4_sol << "\n";
            std::cout << "CVODES: " << cv_sol << "\n";
            std::cout << "CHEB: " << cheb_sol << "\n";
            std::cout << "PS:" << ps_sol(Slice(0, 13)) << "\n";

            //std::cout << "PS:" << ps_sol.size() << "\n";
            //std::cout << "ps_sol: " << ps_sol << "\n";
            /* Reshape to matrix, state at one time equals one column */
            DM ps_sol1 = casadi::DM::reshape(ps_sol, state_size, ps_sol.size1() / state_size);

            DM ps_sol_end = ps_sol1(Slice(0, 13), 0);
            DM ps_sol_0 = ps_sol1(Slice(0, 13), ps_sol1.size2() - 1);
            std::cout << "ps_sol_0: " << ps_sol_0 << "\n";
            std::cout << "ps_sol_end: " << ps_sol_end << "\n";

            sim_states(Slice(1, 1 + state_size), iDatapoint + 1) = ps_sol_end;

        }

        std::cout << "\n"
                     "Seq state_0:   " << sim_states(Slice(0, sim_states.size1()), 0) << "\n";
        std::cout << "Seq state_end: " << sim_states(Slice(0, sim_states.size1()), sim_states.size2() - 1) << "\n";

        std::ofstream trajectory_file(flightDataPath + "states_sim.txt", std::ios::out);
        trajectory_file.precision(15);
        for (int iStep = 0; iStep < sim_states.size2(); ++iStep) {
            //trajectory_file << static_cast<std::vector<double>> (sim_states(Slice(0, sim_states.size1()), i).nonzeros()) << "\n";
            for (int jEntry = 0; jEntry < sim_states.size1(); ++jEntry) {
                trajectory_file << static_cast<double>(sim_states(jEntry, iStep));
                if (jEntry < sim_states.size1() - 1) { trajectory_file << " "; }
            }
            trajectory_file << "\n";
        }
        trajectory_file.close();

        BOOST_CHECK(true);
    }

BOOST_AUTO_TEST_SUITE_END()