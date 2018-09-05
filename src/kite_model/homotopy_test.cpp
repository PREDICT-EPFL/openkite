#include "integrator.h"
#include <fstream>
#include "pseudospectral/chebyshev.hpp"
#include "kiteNMPF.h"

int main(void)
{
    /** Load control signal */
    std::ifstream id_control_file("id_data_control.txt", std::ios::in);
    const int DATA_POINTS = 101;
    const int state_size   = 13;
    const int control_size = 3;

    casadi::DM id_control = casadi::DM::zeros(control_size, DATA_POINTS);

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

    /** create the kite model */
    std::string kite_params_file = "umx_radian.yaml";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_params_file);

    /** no tether */
    kite_props.Tether.Ks = 0.0;
    kite_props.Tether.Kd = 0.0;

    AlgorithmProperties algo_props;
    algo_props.Integrator = CVODES;
    algo_props.sampling_time = 0.02;
    KiteDynamics kite_int(kite_props, algo_props); //integration model
    casadi::Function ode = kite_int.getNumericDynamics();

    /** get dynamics function and state Jacobian */
    casadi::Function DynamicsFunc = kite_int.getNumericDynamics();
    casadi::SX X = kite_int.getSymbolicState();
    casadi::SX U = kite_int.getSymbolicControl();

    /** state bounds */
    casadi::DM INF = casadi::DM::inf(1);
    casadi::DM LBX = casadi::DM::vertcat({2.0, -INF, -INF, -4 * M_PI, -4 * M_PI, -4 * M_PI, -INF, -INF, -INF,
                                          -1.05, -1.05, -1.05, -1.05});
    casadi::DM UBX = casadi::DM::vertcat({INF, INF, INF, 4 * M_PI, 4 * M_PI, 4 * M_PI, INF, INF, INF,
                                          1.05, 1.05, 1.05, 1.05});
    /** control bounds */
    casadi::DM LBU = casadi::DM::vec(id_control);
    casadi::DM UBU = casadi::DM::vec(id_control);

    /** ----------------------------------------------------------------------------------*/
    const int num_segments = 25;
    const int poly_order   = 4;
    const int dimx         = 13;
    const int dimu         = 3;
    const int dimp         = 0;
    const double tf        = 5.0;

    Chebyshev<casadi::SX, poly_order, num_segments, dimx, dimu, dimp> spectral;
    casadi::SX diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, tf);
    diff_constr = diff_constr(casadi::Slice(0, num_segments * poly_order * dimx));

    casadi::SX varx = spectral.VarX();
    casadi::SX varu = spectral.VarU();
    casadi::SX varp = spectral.VarP();

    casadi::SX opt_var = casadi::SX::vertcat(casadi::SXVector{varx, varu, varp});

    casadi::SX lbg = casadi::SX::zeros(diff_constr.size());
    casadi::SX ubg = casadi::SX::zeros(diff_constr.size());

    /** set inequality (box) constraints */
    /** state */
    casadi::SX lbx = casadi::SX::repmat(LBX, num_segments * poly_order + 1, 1);
    casadi::SX ubx = casadi::SX::repmat(UBX, num_segments * poly_order + 1, 1);

    /** control */
    lbx = casadi::SX::vertcat({lbx, LBU});
    ubx = casadi::SX::vertcat({ubx, UBU});

    /** formulate NLP */
    casadi::SXDict NLP;
    casadi::Dict OPTS;
    casadi::DMDict ARG;
    NLP["x"] = opt_var;
    NLP["f"] = 0;
    NLP["g"] = diff_constr;

    OPTS["ipopt.linear_solver"]  = "ma97";
    OPTS["ipopt.print_level"]    = 5;
    OPTS["ipopt.tol"]            = 1e-4;
    OPTS["ipopt.acceptable_tol"] = 1e-4;
    OPTS["ipopt.warm_start_init_point"] = "yes";

    casadi::Function NLP_Solver = casadi::nlpsol("solver", "ipopt", NLP, OPTS);

    /** set default args */
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    casadi::DM feasible_control = (UBU + LBU) / 2;
    casadi::DM init_state = casadi::DM::vertcat({6.1977743e+00,  -2.8407148e-02,   9.1815942e-01,   2.9763089e-01,  -2.2052198e+00,  -1.4827499e-01,
                                         -4.1624807e-01, -2.2601052e+00,   1.2903439e+00,   3.5646195e-02,  -6.9986094e-02,   8.2660637e-01,   5.5727089e-01});
    casadi::DM feasible_state = casadi::DM::repmat(init_state, (num_segments * poly_order + 1), 1);

    ARG["x0"] = casadi::DM::vertcat(casadi::DMVector{feasible_state, feasible_control});

    int idx_in = num_segments * poly_order * dimx;
    int idx_out = idx_in + dimx;
    ARG["lbx"](casadi::Slice(idx_in, idx_out), 0) = init_state;
    ARG["ubx"](casadi::Slice(idx_in, idx_out), 0) = init_state;

    casadi::DMDict res = NLP_Solver(ARG);
    casadi::DM result = res.at("x");

    return 0;
}
