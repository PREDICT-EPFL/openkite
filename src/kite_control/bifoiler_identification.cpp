#include "kiteNMPF.h"
#include "integrator.h"

#define BOOST_TEST_MODULE bifoiler_identification
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include "pseudospectral/chebyshev.hpp"

#include "boat_model.h"

using namespace casadi;

BOOST_AUTO_TEST_SUITE( bifoiler_identification )

BOOST_AUTO_TEST_CASE( bifoiler_id_test )
{
    /** Load identification data */
    std::ifstream id_data_file("id_data_state.txt", std::ios::in);
    std::ifstream id_control_file("id_data_control.txt", std::ios::in);
    const int DATA_POINTS  = 401;
    const int state_size   = 13;
    const int control_size = 4;

    DM id_data    = DM::zeros(state_size, DATA_POINTS);
    DM id_control = DM::zeros(control_size, DATA_POINTS);

    /** load state trajectory */
    if(!id_data_file.fail())
    {
    for(uint i = 0; i < DATA_POINTS; ++i) {
        for(uint j = 0; j < state_size; ++j){
            double entry;
            id_data_file >> entry;
            id_data(j,i) = entry;
        }
    }
    }
    else
    {
        std::cout << "Could not open : id state data file \n";
        id_data_file.clear();
    }

    /** load control data */
    if(!id_control_file.fail())
    {
    for(uint i = 0; i < DATA_POINTS; ++i){
        for(uint j = 0; j < control_size; ++j){
            double entry;
            id_control_file >> entry;
            /** put in reverse order to comply with Chebyshev method */
            id_control(j,DATA_POINTS - 1 - i) = entry;
        }
    }
    }
    else
    {
        std::cout << "Could not open : id control data file \n";
        id_control_file.clear();
    }


    /** define boat dynamics */
    std::string boat_config_file = "config.yaml";
    BoatProperties boat_props = BoatProperties::Load(boat_config_file);

    bifoiler::BoatDynamics boat(boat_props, true);
    bifoiler::BoatDynamics boat_int(boat_props);    //integration model
    Function ode = boat_int.getNumericDynamics();

    /** get dynamics function and state Jacobian */
    Function DynamicsFunc = boat.getNumericDynamics();
    SX X = boat.getSymbolicState();
    SX U = boat.getSymbolicControl();
    SX P = boat.getSymbolicParams();

    std::cout << "OK \n";


    /** state bounds */
    DM LBX = DM::vertcat({0.5, -DM::inf(1), -DM::inf(1), -4 * M_PI, -4 * M_PI, -4 * M_PI, -DM::inf(1), -DM::inf(1), -DM::inf(1),
                          -1.05, -1.05, -1.05, -1.05});
    DM UBX = DM::vertcat({DM::inf(1), DM::inf(1), DM::inf(1), 4 * M_PI, 4 * M_PI, 4 * M_PI, DM::inf(1), DM::inf(1), DM::inf(1),
                          1.05, 1.05, 1.05, 1.05});
    /** control bounds */
    DM LBU = DM::vec(id_control);
    DM UBU = DM::vec(id_control);

    /** parameter bounds */
    YAML::Node config = YAML::LoadFile("config.yaml");
    double CL0     = config["hydrodynamic"]["CL0"].as<double>();
    double CLa_tot = config["hydrodynamic"]["CLa_total"].as<double>();
    double CD0_tot = config["hydrodynamic"]["CD0_total"].as<double>();

    double CYb = config["hydrodynamic"]["CYb"].as<double>();
    double Cm0 = config["hydrodynamic"]["Cm0"].as<double>();
    double Cma = config["hydrodynamic"]["Cma"].as<double>();
    double Cnb = config["hydrodynamic"]["Cnb"].as<double>();
    double Clb = config["hydrodynamic"]["Clb"].as<double>();

    double CLq = config["hydrodynamic"]["CLq"].as<double>();
    double Cmq = config["hydrodynamic"]["Cmq"].as<double>();
    double CYr = config["hydrodynamic"]["CYr"].as<double>();
    double Cnr = config["hydrodynamic"]["Cnr"].as<double>();
    double Clr = config["hydrodynamic"]["Clr"].as<double>();

    double CYp = config["hydrodynamic"]["CYp"].as<double>();
    double Clp = config["hydrodynamic"]["Clp"].as<double>();
    double Cnp = config["hydrodynamic"]["Cnp"].as<double>();

    double CXdf = config["hydrodynamic"]["CXdf"].as<double>();
    double CYdr = config["hydrodynamic"]["CYdr"].as<double>();
    double CZdf = config["hydrodynamic"]["CZdf"].as<double>();
    double Clda = config["hydrodynamic"]["CLda"].as<double>();
    double Cldr = config["hydrodynamic"]["CLdr"].as<double>();
    double Cmdf = config["hydrodynamic"]["CMdf"].as<double>();
    double Cnda = config["hydrodynamic"]["CNda"].as<double>();
    double Cndr = config["hydrodynamic"]["CNdr"].as<double>();


    DM REF_P = DM::vertcat({CL0, CLa_tot, CD0_tot, CYb, Cm0, Cma, Cnb, Clb, CLq, Cmq, CYr, Cnr,
                            Clr, CYp, Clp, Cnp, CXdf, CYdr, CZdf, Clda, Cldr, Cmdf, Cnda, Cndr});
    DM LBP = REF_P; DM UBP = REF_P;
    LBP = -DM::inf(24);
    UBP = DM::inf(24);


    /** ----------------------------------------------------------------------------------*/
    const int num_segments = 20;
    const int poly_order   = 20;
    const int dimx         = 13;
    const int dimu         = 4;
    const int dimp         = 24;
    const double tf        = 19.55;

    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> spectral;
    SX diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, tf);

    diff_constr = diff_constr(casadi::Slice(0, num_segments * poly_order * dimx));

    SX varx = spectral.VarX();
    SX varu = spectral.VarU();
    SX varp = spectral.VarP();

    SX opt_var = SX::vertcat(SXVector{varx, varu, varp});

    SX lbg = SX::zeros(diff_constr.size());
    SX ubg = SX::zeros(diff_constr.size());

    /** set inequality (box) constraints */
    /** state */
    SX lbx = SX::repmat(LBX, num_segments * poly_order + 1, 1);
    SX ubx = SX::repmat(UBX, num_segments * poly_order + 1, 1);

    /** control */
    lbx = SX::vertcat({lbx, LBU});
    ubx = SX::vertcat({ubx, UBU});

    /** parameters */
    lbx = SX::vertcat({lbx, LBP});
    ubx = SX::vertcat({ubx, UBP});


    DM Q  = SX::diag(SX({1e3, 1e3, 1e3,  1e3, 1e3, 1e3,  1e1, 1e1, 1e1,  1e1, 1e1, 1e1, 1e1})); //good one as well
    //DM Q = 1e1 * DM::eye(13);
    double alpha = 0.0;


    SX fitting_error = 0;
    SX varx_ = SX::reshape(varx, 13, DATA_POINTS);
    for (uint j = 0; j < DATA_POINTS; ++j)
    {
        SX measurement = id_data(Slice(0, id_data.size1()), j);
        SX error = measurement - varx_(Slice(0, varx_.size1()), varx_.size2() - j - 1);
        fitting_error = fitting_error + (1 / DATA_POINTS) * SX::sumRows( SX::mtimes(Q, pow(error, 2)) );
    }

    /** add regularisation */
    fitting_error = fitting_error + alpha * SX::dot(varp - SX({REF_P}), varp - SX({REF_P}));

    /** alternative approximation */
    SX x = SX::sym("x",13);
    SX y = SX::sym("y",13);
    SX cost_function = SX::sumRows( SX::mtimes(Q, pow(x - y, 2)) );
    Function IdCost = Function("IdCost",{x,y}, {cost_function});
    SX fitting_error2 = spectral.CollocateIdCost(IdCost, id_data, 0, tf);
    fitting_error2 = fitting_error2 + alpha * SX::dot(varp - SX({REF_P}), varp - SX({REF_P}));

    /** formulate NLP */
    SXDict NLP;
    Dict OPTS;
    DMDict ARG;
    NLP["x"] = opt_var;
    NLP["f"] = fitting_error;
    NLP["g"] = diff_constr;

    OPTS["ipopt.linear_solver"]  = "ma97";
    OPTS["ipopt.print_level"]    = 5;
    OPTS["ipopt.tol"]            = 1e-4;
    OPTS["ipopt.acceptable_tol"] = 1e-4;
    OPTS["ipopt.warm_start_init_point"] = "yes";
    //OPTS["ipopt.max_iter"]       = 20;

    Function NLP_Solver = nlpsol("solver", "ipopt", NLP, OPTS);

    /** set default args */
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    /** provide an initial guess from the integrator */
    casadi::DMDict props;
    props["scale"] = 0;
    props["P"] = casadi::DM::diag(casadi::DM({0.1, 1/3.0, 1/3.0, 1/2.0, 1/5.0, 1/2.0, 1/3.0, 1/3.0, 1/3.0, 1.0, 1.0, 1.0, 1.0}));
    props["R"] = casadi::DM::diag(casadi::DM({1.0, 1/0.15, 1/0.2618, 1/0.2618}));

    PSODESolver<poly_order,num_segments,dimx,dimu>ps_solver(ode, tf, props);

    DM x0 = id_data(Slice(0, id_data.size1()), 0);
    //std::cout << DM::vec(id_data) << "\n";
    //DM x0 = DM::repmat(id_data(Slice(0, 13), (NumSegments * PolyOrder + 1), 1);
    DMDict solution = ps_solver.solve_trajectory(x0, LBU, true);
    DM feasible_state = solution.at("x");

    /**
    DM feasible_state = DM::reshape(id_data, dimx * (num_segments * poly_order + 1), 1);
    //DM feasible_state = DM::repmat(id_data(Slice(0, id_data.size1()), 0), (num_segments * poly_order + 1), 1);
    DM feasible_control = (UBU + LBU) / 2;

    ARG["x0"] = DM::vertcat(DMVector{feasible_state, feasible_control, REF_P});
    //ARG["x0"] = DM::vertcat(DMVector{feasible_state, REF_P});
    //ARG["lam_g0"] = solution.at("lam_g");
    //ARG["lam_x0"] = DM::vertcat({solution.at("lam_x"), DM::zeros(REF_P.size1())});

    int idx_in = num_segments * poly_order * dimx;
    int idx_out = idx_in + dimx;
    ARG["lbx"](Slice(idx_in, idx_out), 0) = x0;
    ARG["ubx"](Slice(idx_in, idx_out), 0) = x0;

    DMDict res = NLP_Solver(ARG);
    DM result = res.at("x");

    DM new_params = result(Slice(result.size1() - varp.size1(), result.size1()));
    std::vector<double> new_params_vec = new_params.nonzeros();

    DM trajectory = result(Slice(0, varx.size1()));
    //DM trajectory = DM::reshape(traj, DATA_POINTS, dimx );
    std::ofstream trajectory_file("estimated_trajectory.txt", std::ios::out);

    if(!trajectory_file.fail())
    {
        for (int i = 0; i < trajectory.size1(); i = i + dimx)
        {
            std::vector<double> tmp = trajectory(Slice(i, i + dimx),0).nonzeros();
            for (uint j = 0; j < tmp.size(); j++)
            {
                trajectory_file << tmp[j] << " ";
            }
            trajectory_file << "\n";
        }
    }
    trajectory_file.close();
    */

    /** update parameter file */
    /**
    config["hydrodynamic"]["CL0"] = new_params_vec[0];
    config["hydrodynamic"]["CLa_total"] = new_params_vec[1];
    config["hydrodynamic"]["CD0_total"] = new_params_vec[2];

    config["hydrodynamic"]["CYb"] = new_params_vec[3];
    config["hydrodynamic"]["Cm0"] = new_params_vec[4];
    config["hydrodynamic"]["Cma"] = new_params_vec[5];
    config["hydrodynamic"]["Cnb"] = new_params_vec[6];
    config["hydrodynamic"]["Clb"] = new_params_vec[7];

    config["hydrodynamic"]["CLq"] = new_params_vec[8];
    config["hydrodynamic"]["Cmq"] = new_params_vec[9];
    config["hydrodynamic"]["CYr"] = new_params_vec[10];
    config["hydrodynamic"]["Cnr"] = new_params_vec[11];
    config["hydrodynamic"]["Clr"] = new_params_vec[12];

    config["hydrodynamic"]["CYp"] = new_params_vec[13];
    config["hydrodynamic"]["Clp"] = new_params_vec[14];
    config["hydrodynamic"]["Cnp"] = new_params_vec[15];

    config["hydrodynamic"]["CXdf"] = new_params_vec[16];
    config["hydrodynamic"]["CYdr"] = new_params_vec[17];
    config["hydrodynamic"]["CZdf"] = new_params_vec[18];
    config["hydrodynamic"]["CLda"] = new_params_vec[19];
    config["hydrodynamic"]["CLdr"] = new_params_vec[20];
    config["hydrodynamic"]["CMdf"] = new_params_vec[21];
    config["hydrodynamic"]["CNda"] = new_params_vec[22];
    config["hydrodynamic"]["CNdr"] = new_params_vec[23];

    std::ofstream fout("config_id.yaml");
    fout << config;
    */

    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_SUITE_END()
