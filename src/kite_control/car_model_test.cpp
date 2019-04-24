#include "mobile_robot.h"
#include "kinematic_robot.hpp"
#include "integrator.h"
#include "kite.h"
#include "iomanip"

using namespace casadi;


int main(void)
{
    //KinematicRobotProperties robo_props(2.9);
    //KinematicRobot robot(robo_props);
    KinematicRobot robot;
    Function ode = robot.getDynamics();

    int dimx = 6;

    /** compare three ode solvers */
    Dict opts;
    opts["tf"]         = 0.01;
    opts["poly_order"] = 5;
    opts["tol"]        = 1e-3;
    opts["method"] = IntType::CVODES;

    double tf = 0.01;
    ODESolver cvodes_solver(ode, opts);
    casadi::DMDict props;

    /** scaling */
    DM P = 0.1 * DM::eye(8);
    DM R = DM::diag(DM({1,1/4000.0, 1/4000.0}));
    props["scale"] = 0;
    props["P"] = P;
    props["R"] = R;
    props["reg_lag"] = 0;
    PSODESolver<4,2,6,3>ps_solver(ode, 2.0, props);

    /** read the telemetry log */
    DMVector telemetry;
    std::ifstream telemetry_file("mpc_log2.txt", std::ios::in);

    /** load control data */
    if(!telemetry_file.fail())
    {
        while(!telemetry_file.eof())
        {
            std::vector<double> line;
            double entry;
            for (int i = 0; i < 26; ++i)
            {
                telemetry_file >> entry;
                line.push_back(entry);
            }
            telemetry.push_back(DM({line}).T());
        }
    }
    else
    {
        std::cout << "Could not open : telemetry file \n";
        telemetry_file.clear();
    }

    DM telemetry_mat = DM::vertcat(telemetry);
    DM control_vec = telemetry_mat(Slice(0, telemetry_mat.size1()), Slice(1, 4));

    /** solve a problem */
    DM cv_sol, force, ode_val;
    DM init_state = telemetry[0](Slice(4, 4 + dimx));
    //DM control    = control_vec(0,Slice(0,3));
    //DM control = DM::zeros(2);
    //control[0] = 5; // control_log[0];
    //control[1] = 0.2; //control_log[2] - control_log[1];
    //init_state[0] = 0.5;
    //control[2] = 4000;
    //control = DM({5, 0.2});
    //control  =  DM({0.1, 0, 4000});
    DM control = DM({0.0, 0.0, 100});

    DM control_traj = DM::repmat(control, 9, 1);
    //init_state = DM({0,0,0, 0,0,0,0.0,0.0}).T();

    std::chrono::time_point<std::chrono::system_clock> start = kite_utils::get_time();
    ode_val = ode(DMVector{init_state, control})[0];
    DMDict ps_sol   = ps_solver.solve_trajectory(init_state.T(), control_traj);
    cv_sol = cvodes_solver.solve(init_state, control, tf);
    force = robot.TraceFunction(DMVector{init_state, control})[0];
    std::chrono::time_point<std::chrono::system_clock> cv_stop = kite_utils::get_time();

    auto cv_duration = std::chrono::duration_cast<std::chrono::microseconds>(cv_stop - start);

    std::cout << "CVODES solve time: " << std::setprecision(6)
              << static_cast<double>(cv_duration.count()) * 1e-6 << " [seconds]" << "\n";

    std::cout << "CVODES: " <<  cv_sol << "\n";
    std::cout << "ODE value: " << ode_val << "\n";
    std::cout << "Forces: " << force << "\n";
    std::cout << "X0: " << init_state << "\n";
    std::cout << "U: " << control << "\n";

    /** simulation loop */
    DMVector sim_log;
    DM time = 0;
    DM two_zeros = DM::zeros(2);
    for(int i = 0; i < 0; ++i)
    {
        control = control_vec(i,Slice(0,3));
        //DM fake_state = telemetry[i](Slice(4, 4+8));

        std::cout << "Iteration: " << i << "\n";
        cv_sol = cvodes_solver.solve(init_state, control, tf);
        DM dynamics = ode(DMVector{cv_sol, control})[0];
        DM forces = robot.TraceFunction(DMVector{init_state, control})[0];

        init_state = cv_sol;
        time += tf;
        //DM log_msg = DM::vertcat(DMVector{time, control.T(), cv_sol, dynamics, forces});
        //simple dynamics version
        DM log_msg = DM::vertcat(DMVector{time, control.T(), cv_sol, two_zeros, dynamics, two_zeros, forces});
        sim_log.push_back(log_msg.T());
    }


    std::cout << "Simulation is over \n";

    /** save simulation to a file */
    std::ofstream simulation_file("mpc_simulation.txt", std::ios::out);
    if(!simulation_file.fail())
    {
        for (int i = 0; i < sim_log.size(); ++i)
        {
            simulation_file << sim_log[i].nonzeros() << "\n";
        }
    }
    simulation_file.close();

    /** static model test [accelerations] */
    DM static_init = DM({8.3086,   -0.3783,   -0.4190,  -81.6415,  179.5490,   -2.0047});
    DM static_input = DM({-0.140761, 0, 2322.14});
    DM acc = ode(DMVector{static_init, static_input})[0];

    std::cout << "deriv: " << acc.nonzeros() << "\n";

    return 0;
}
