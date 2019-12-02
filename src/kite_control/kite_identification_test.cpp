#include "kiteNMPF.h"
#include "integrator.h"
#include <fstream>
#include "pseudospectral/chebyshev.hpp"

using namespace casadi;


void set_parameter_bounds(DM &LBP, DM &UBP, const int &index,
                          const DM &REF_P, const double &LB_factor, const double &UB_factor) {

    LBP(index) = REF_P(index) + LB_factor * fabs(REF_P(index));
    UBP(index) = REF_P(index) + UB_factor * fabs(REF_P(index));

}

void set_parameter_bounds(DM &LBP, DM &UBP, const int &index,
                          const double &absoluteLB, const double &absoluteUB) {

    LBP(index) = absoluteLB;
    UBP(index) = absoluteUB;

}

int main(void) {

    std::string flightDataDir = "/home/johannes/identification/processed_flight_data/";

    struct FlightMetaData {
        /// 1. ///
        int session = 2;
        int number = 1;
        int seq = 4;

        /// 2. ///
        //std::string maneuver = "identRoll";
        std::string maneuver = "identPitch";

        /// 3. ///
        std::string resampleMethod = "cheb";
        //std::string resampleMethod = "equal";
    } flight;

    /// 4. ///
    const int dimp = 1 + 9; // 1 general parameter + // lon: 9, lat: 16 identification parameters

    /// 5. ///
    const int DATA_POINTS = 127;
    const int num_segments = 42;
    const int poly_order = 3;

    /** COMMENT / UNCOMMENT THE SUITING PARAMETERS IN KITE.CPP ! */

    std::string flightDataPath = flightDataDir
                                 + "session_" + std::to_string(flight.session)
                                 + "/flight_" + std::to_string(flight.number);
    if (!flight.maneuver.empty()) flightDataPath.append("/" + flight.maneuver);
    flightDataPath.append("/seq_" + std::to_string(flight.seq) + "/");

    std::cout << flightDataPath << "\n";

    /** define kite dynamics */
    std::string kite_params_file = "/home/johannes/identification/easy_glider_4.yaml";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_params_file);

    /** Load wind data **/
    std::ifstream id_wind_file(flightDataPath + "wind.txt", std::ios::in);
    if (!id_wind_file.fail()) {
        int windFrom_deg;
        double windSpeed;
        id_wind_file >> windFrom_deg;
        id_wind_file >> windSpeed;

        kite_props.Wind.WindFrom_deg = windFrom_deg;
        kite_props.Wind.WindSpeed = windSpeed;
    } else {
        std::cout << "Could not open : id wind data file \n";
        id_wind_file.clear();
    }

    /** Load identification data */
    std::ifstream id_data_file(flightDataPath + flight.resampleMethod + "_states.txt", std::ios::in);
    std::ifstream id_control_file(flightDataPath + flight.resampleMethod + "_controls.txt", std::ios::in);
    const int dimx = 13;  // v(3) w(3) r(3) q(4)
    const int dimu = 4; // T elev rud ail

    DM id_data = DM::zeros(1 + dimx, DATA_POINTS);      // Time + States
    DM id_control = DM::zeros(1 + dimu, DATA_POINTS);   // Time + Controls

    /** load state trajectory */
    if (!id_data_file.fail()) {
        for (uint i = 0; i < DATA_POINTS; ++i) {
            for (uint j = 0; j < (1 + dimx); ++j) {       // Time + States
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
            for (uint j = 0; j < (1 + dimu); ++j) {     // Time + Controls
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
    algo_props.Integrator = CVODES;
    algo_props.sampling_time = 0.02;
    //algo_props.sampling_time = static_cast<double>(id_data(0, 1) - id_data(0, 0));

    KiteDynamics kite(kite_props, algo_props, true);
    KiteDynamics kite_int(kite_props, algo_props); //integration model
    Function ode = kite_int.getNumericDynamics();

    /** get dynamics function and state Jacobian */
    Function DynamicsFunc = kite.getNumericDynamics();
    SX X = kite.getSymbolicState();
    SX U = kite.getSymbolicControl();
    SX P = kite.getSymbolicParameters();

    /** state bounds */
    DM LBX = DM::vertcat({2.0, -50, -50,
                          -4 * M_PI, -4 * M_PI, -4 * M_PI,
                          -300, -300, -300,
                          -1.05, -1.05, -1.05, -1.05});     // for infinity use -DM::inf(1) and DM::inf(1)
    DM UBX = DM::vertcat({50, 50, 50,
                          4 * M_PI, 4 * M_PI, 4 * M_PI,
                          300, 300, 300,
//                          10, 10, 10,
                          1.05, 1.05, 1.05, 1.05});

    /** control bounds */
    DM LBU = DM::vec(
            id_control(Slice(1, id_control.size1()), Slice(0, id_control.size2()))); // without Time (first row)
    DM UBU = DM::vec(
            id_control(Slice(1, id_control.size1()), Slice(0, id_control.size2()))); // without Time (first row)

    /** parameter bounds */
    DM REF_P;
    DM LBP;
    DM UBP;

    /** Cost matrix **/
    DM Q;

    /** parameter bounds preparation */
    /** --------------------- **/
    /** Wind properties       **/
    /** --------------------- **/
    double windFrom_deg = kite_props.Wind.WindFrom_deg;
    double windSpeed = kite_props.Wind.WindSpeed;

    /** --------------------- **/
    /** Geometric parameters  **/
    /** --------------------- **/
    double imuPitchOffset = kite_props.Geometry.ImuPitchOffset;

    REF_P = DM::vertcat({REF_P, imuPitchOffset}); // 1 general parameter

    //double b = kite_props.Geometry.WingSpan;
    //double c = kite_props.Geometry.MAC;
    //double AR = kite_props.Geometry.AspectRatio;
    //double S = kite_props.Geometry.WingSurfaceArea;
    //double lam = kite_props.Geometry.TaperRatio;
    //double St = kite_props.Geometry.HTailsurface;
    //double lt = kite_props.Geometry.TailLeverArm;
    //double Sf = kite_props.Geometry.FinSurfaceArea;
    //double lf = kite_props.Geometry.FinLeverArm;
    //double Xac = kite_props.Geometry.AerodynamicCenter;
    /** @todo: get rid of all magic numbers **/
    //double Xcg = 0.031/c;               /** Center of Gravity (CoG) wrt leading edge [1/c] **/
    //double Vf = (Sf * lf) / (S * b);    /** fin volume coefficient []                      **/
    //double Vh = (lt * St) / (S * c);    /** horizontal tail volume coefficient []          **/

    /** --------------------------- **/
    /** Mass and inertia parameters **/
    /** --------------------------- **/
    //double Mass = kite_props.Inertia.Mass;
    //double Ixx = kite_props.Inertia.Ixx;
    //double Iyy = kite_props.Inertia.Iyy;
    //double Izz = kite_props.Inertia.Izz;
    //double Ixz = kite_props.Inertia.Ixz;

    /** ------------------------------- **/
    /** Static aerodynamic coefficients **/
    /** ------------------------------- **/
    double CL0 = kite_props.Aerodynamics.CL0;
    //double CL0_t = kite_props.Aerodynamics.CL0_tail;
    double CLa_tot = kite_props.Aerodynamics.CLa_total;
    //double CLa_w = kite_props.Aerodynamics.CLa_wing;
    //double CLa_t = kite_props.Aerodynamics.CLa_tail;
    //double e_o = kite_props.Aerodynamics.e_oswald;

    double CD0_tot = kite_props.Aerodynamics.CD0_total;
    //double CD0_w = kite_props.Aerodynamics.CD0_wing;
    //double CD0_t = kite_props.Aerodynamics.CD0_tail;
    double CYb = kite_props.Aerodynamics.CYb;
    //double CYb_vt = kite_props.Aerodynamics.CYb_vtail;
    double Cm0 = kite_props.Aerodynamics.Cm0;
    double Cma = kite_props.Aerodynamics.Cma;
    double Cn0 = kite_props.Aerodynamics.Cn0;
    double Cnb = kite_props.Aerodynamics.Cnb;
    double Cl0 = kite_props.Aerodynamics.Cl0;
    double Clb = kite_props.Aerodynamics.Clb;
    //double dw = CLa_tot / (pi * e_o * AR);             /** downwash acting at the tail [] **/

    double CLq = kite_props.Aerodynamics.CLq;
    double Cmq = kite_props.Aerodynamics.Cmq;
    double CYr = kite_props.Aerodynamics.CYr;
    double Cnr = kite_props.Aerodynamics.Cnr;
    double Clr = kite_props.Aerodynamics.Clr;
    double CYp = kite_props.Aerodynamics.CYp;
    double Clp = kite_props.Aerodynamics.Clp;
    double Cnp = kite_props.Aerodynamics.Cnp;

    /** ------------------------------ **/
    /** Aerodynamic effects of control **/
    /** ------------------------------ **/
    double CLde = kite_props.Aerodynamics.CLde;
    double CYdr = kite_props.Aerodynamics.CYdr;
    double Cmde = kite_props.Aerodynamics.Cmde;
    double Cndr = kite_props.Aerodynamics.Cndr;
    double Cldr = kite_props.Aerodynamics.Cldr;
    //double CDde = kite_props.Aerodynamics.CDde;
    double Clda = kite_props.Aerodynamics.Clda;
    double Cnda = kite_props.Aerodynamics.Cnda;

    //double CL_daoa = -2 * CLa_t * Vh * dw;
    //double Cm_daoa = -2 * CLa_t * Vh * (lt/c) * dw;

    /** ------------------------------ **/
    /**        Tether parameters       **/
    /** ------------------------------ **/
    //    double Lt = kite_props.Tether.length;
    //    double Ks = kite_props.Tether.Ks;
    //    double Kd = kite_props.Tether.Kd;
    //    double rx = kite_props.Tether.rx;
    //    double ry = kite_props.Tether.ry;
    //    double rz = kite_props.Tether.rz;

    // LBP(21] = 2.65;    UBP(21] = 2.75;   // tether length
    // LBP(22] = 150.0;  UBP(22] = 150.0;  // Ks
    // LBP(23] = 0.0;    UBP(23] = 10;   // Kd
    // LBP(24] = 0.0;    UBP(24] = 0.0;   // rx
    // LBP(25] = 0.0;    UBP(25] = 0.0;  // rz

    //DM REF_P = DM::vertcat({Lt, Ks, Kd, rx, ry, rz});

    /** LONGITUDINAL IDENTIFICATION PARAMETERS ---------------------------------------------------------------------- */
    if (flight.maneuver == "identPitch") {
        int b = REF_P.size1(); // parameter size before
        REF_P = DM::vertcat({REF_P,
                             CL0, CLa_tot,
                             CD0_tot, Cm0, Cma,
                             CLq, Cmq,
                             CLde, Cmde,
                            }); // 9 longitudinal parameters

        /** parameter bounds */
        LBP = REF_P;
        UBP = REF_P;

        set_parameter_bounds(LBP, UBP, b + 0, REF_P, -0.1, 0.1); // CL0
        set_parameter_bounds(LBP, UBP, b + 1, REF_P, -0.05, 0.1); // CLa_tot

        set_parameter_bounds(LBP, UBP, b + 2, REF_P, -0.1, 0.25); // CD0_tot
        set_parameter_bounds(LBP, UBP, b + 3, REF_P, -0.5, 0.5); // Cm0
        set_parameter_bounds(LBP, UBP, b + 4, REF_P, -0.1, 0.3); // Cma

        set_parameter_bounds(LBP, UBP, b + 5, REF_P, -0.2, 0.2); // CLq
        set_parameter_bounds(LBP, UBP, b + 6, REF_P, -0.3, 0.3); // Cmq

        set_parameter_bounds(LBP, UBP, b + 7, REF_P, -0.5, 0.5); // CLde
        set_parameter_bounds(LBP, UBP, b + 8, REF_P, -0.5, 0.5); // Cmde

        Q = SX::diag(SX({1e2, 0, 1e2,
                         0, 1e2, 0,
                         1e1, 1e1, 1e2,
                         1e2, 1e2, 1e2, 1e2}));
    }
    /** END OF LONGITUDINAL IDENTIFICATION PARAMETERS --------------------------------------------------------------- */

    /** LATERAL IDENTIFICATION PARAMETERS --------------------------------------------------------------------------- */
    if (flight.maneuver == "identRoll") {
        int b = REF_P.size1(); // parameter size before
        REF_P = DM::vertcat({REF_P,
                             CYb, Cn0, Cnb, Cl0, Clb,
                             CYr, Cnr, Clr, CYp, Clp, Cnp,
                             CYdr, Cndr, Cldr, Clda, Cnda,
                            }); // 16 lateral parameters

        /** parameter bounds */
        LBP = REF_P;
        UBP = REF_P;

        set_parameter_bounds(LBP, UBP, b + 0, REF_P, -0.5, 0.5);  // CYb
        set_parameter_bounds(LBP, UBP, b + 1, REF_P, -0.5, 0.5);  // Cn0
        set_parameter_bounds(LBP, UBP, b + 2, REF_P, -0.5, 0.5);  // Cnb
        set_parameter_bounds(LBP, UBP, b + 3, REF_P, -0.5, 0.5);  // Cl0
        set_parameter_bounds(LBP, UBP, b + 4, REF_P, -0.5, 0.5);  // Clb

        set_parameter_bounds(LBP, UBP, b + 5, REF_P, -0.3, 0.3);  // CYr
        set_parameter_bounds(LBP, UBP, b + 6, REF_P, -0.5, 0.5);  // Cnr
        set_parameter_bounds(LBP, UBP, b + 7, REF_P, -0.5, 0.5);  // Clr
        set_parameter_bounds(LBP, UBP, b + 8, REF_P, -0.5, 0.5);  // CYp
        set_parameter_bounds(LBP, UBP, b + 9, REF_P, -0.5, 0.5);  // Clp
        set_parameter_bounds(LBP, UBP, b + 10, REF_P, -0.3, 1.0);  // Cnp

        set_parameter_bounds(LBP, UBP, b + 11, REF_P, -0.5, 0.5);  // CYdr
        set_parameter_bounds(LBP, UBP, b + 12, REF_P, -0.5, 0.5);  // Cndr
        set_parameter_bounds(LBP, UBP, b + 13, REF_P, -0.5, 0.5);  // Cldr
        set_parameter_bounds(LBP, UBP, b + 14, REF_P, -0.5, 0.5);  // Clda
        set_parameter_bounds(LBP, UBP, b + 15, REF_P, -0.5, 0.5);  // Cnda

        Q = SX::diag(SX({0, 1e2, 0,
                         1e2, 0, 1e2,
                         1e1, 1e1, 1e2,
                         1e2, 1e2, 1e2, 1e2}));
    }
    /** END OF LATERAL IDENTIFICATION PARAMETERS -------------------------------------------------------------------- */

    /** General parameter bounds */
    set_parameter_bounds(LBP, UBP, 0, -7.0 * M_PI / 180.0, 7.0 * M_PI / 180.0);  // imuOffsetPitch

    std::cout << "OK: Parameter preparation\n";

    /** ----------------------------------------------------------------------------------*/
//    const int num_segments = 42;
//    const int poly_order = 3;
//    const int dimx = 13;
//    const int dimu = 4;
    //const double tf = 10.0;
    //int dimx = static_cast<int>(REF_P.size1());
    //int dimu = static_cast<int>(U.size1());
    //int dimp = static_cast<int>(REF_P.size1());
    double tf = static_cast<double>(id_data(0, id_data.size2() - 1) - id_data(0, 0));

    std::cout << "tf = " << tf << "\n";

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

//    DM Q = SX::diag(SX({1e3, 1e2, 1e2,
//                        1e2, 1e2, 1e2,
//                        1e1, 1e1, 1e2,
//                        1e2, 1e2, 1e2, 1e2}));

    //DM Q = 1e1 * DM::eye(13);
    double alpha = 100.0;

    SX fitting_error = 0;
    SX varx_ = SX::reshape(varx, dimx, DATA_POINTS);

    for (uint j = 0; j < DATA_POINTS; ++j) {
        SX measurement = id_data(Slice(1, id_data.size1()), j);
        SX error = measurement - varx_(Slice(0, varx_.size1()), varx_.size2() - j - 1);
        fitting_error += static_cast<double>(1.0 / DATA_POINTS) * SX::sum1(SX::mtimes(Q, pow(error, 2)));
    }

    /** add regularisation */
    // fitting_error = fitting_error + alpha * SX::dot(varp - SX({REF_P}), varp - SX({REF_P}));

    /** alternative approximation */
    /**
    SX x = SX::sym("x", dimx);
    SX y = SX::sym("y", dimx);
    SX cost_function = SX::sumRows( SX::mtimes(Q, pow(x - y, 2)) );
    Function IdCost = Function("IdCost",{x,y}, {cost_function});
    SX fitting_error2 = spectral.CollocateIdCost(IdCost, id_data, 0, tf);
    fitting_error2 = fitting_error2 + alpha * SX::dot(varp - SX({REF_P}), varp - SX({REF_P}));
    */

    /** formulate NLP */
    SXDict NLP;
    Dict OPTS;
    DMDict ARG;
    NLP["x"] = opt_var;
    NLP["f"] = fitting_error;
    NLP["g"] = diff_constr;

    OPTS["ipopt.linear_solver"] = "mumps";
    OPTS["ipopt.print_level"] = 5;
    OPTS["ipopt.tol"] = 1e-4;
    OPTS["ipopt.acceptable_tol"] = 1e-4;
    OPTS["ipopt.warm_start_init_point"] = "yes";
    //OPTS["ipopt.max_iter"]       = 20;

    Function NLP_Solver = nlpsol("solver", "ipopt", NLP, OPTS);

    std::cout << "OK: Solver set up\n";

    /** set default args */
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    DMDict solution;
    DM feasible_state;
    DM init_state = id_data(Slice(1, id_data.size1()), 0);

    /** if the solutions available load them from file */
    if (kite_utils::file_exists(flightDataPath + "id_x0.txt")) {
        DM sol_x = kite_utils::read_from_file(flightDataPath + "id_x0.txt");
        ARG["x0"] = DM::vertcat(DMVector{sol_x, REF_P});
        feasible_state = sol_x;

        std::cout << "Initial guess loaded from a file \n";
    } else {
        /** otherwise, provide initial guess from integrator */
        casadi::DMDict props;
        props["scale"] = 0;
        props["P"] = casadi::DM::diag(casadi::DM(
                {0.1, 1 / 3.0, 1 / 3.0,
                 1 / 2.0, 1 / 5.0, 1 / 2.0,
                 1 / 3.0, 1 / 3.0, 1 / 3.0,
                 1.0, 1.0, 1.0, 1.0}));
        PSODESolver<poly_order, num_segments, dimx, dimu> ps_solver(ode, tf, props);

        DM init_control = DM({0.1, 0.0, 0.0, 0.0});
        init_control = casadi::DM::repmat(init_control, (num_segments * poly_order + 1), 1);
        std::cout << "init_state: " << init_state << " init_control: " << init_control << "\n";
        solution = ps_solver.solve_trajectory(init_state, init_control, true);
        feasible_state = solution.at("x");
        ARG["x0"] = casadi::DM::vertcat(casadi::DMVector{feasible_state, REF_P});
    }

    if (kite_utils::file_exists(flightDataPath + "id_lam_g.txt")) {
        ARG["lam_g0"] = kite_utils::read_from_file(flightDataPath + "id_lam_g.txt");
    } else {
        ARG["lam_g0"] = solution.at("lam_g");
    }

    if (kite_utils::file_exists(flightDataPath + "id_lam_x.txt")) {
        DM sol_lam_x = kite_utils::read_from_file(flightDataPath + "id_lam_x.txt");
        ARG["lam_x0"] = DM::vertcat({sol_lam_x, DM::zeros(REF_P.size1())});
    } else {
        ARG["lam_x0"] = DM::vertcat({solution.at("lam_x"), DM::zeros(REF_P.size1())});
    }

//    /** write initial trajectory to a file */
//    std::ofstream trajectory_file(flightDataPath + "integrated_trajectory.txt", std::ios::out);
//    if (!trajectory_file.fail()) {
//        for (int i = 0; i < varx.size1(); i = i + 13) {
//            std::vector<double> tmp = feasible_state(Slice(i, i + 13), 0).nonzeros();
//            for (uint j = 0; j < tmp.size(); j++) {
//                trajectory_file << tmp[j] << " ";
//            }
//            trajectory_file << "\n";
//        }
//    }
//    trajectory_file.close();

    int idx_in = num_segments * poly_order * dimx;
    int idx_out = idx_in + dimx;
    ARG["lbx"](Slice(idx_in, idx_out), 0) = init_state;
    ARG["ubx"](Slice(idx_in, idx_out), 0) = init_state;

    /** solve the identification problem */
    DMDict res = NLP_Solver(ARG);
    DM result = res.at("x");
    DM lam_x = res.at("lam_x");

    DM new_params = result(Slice(result.size1() - varp.size1(), result.size1()));
    DM param_sens = lam_x(Slice(lam_x.size1() - varp.size1(), lam_x.size1()));

    std::cout << "PARAMETER SENSITIVITIES: " << param_sens << "\n";

    std::vector<double> new_params_vec = new_params.nonzeros();

//    DM trajectory = result(Slice(0, varx.size1()));
//    //DM trajectory = DM::reshape(traj, DATA_POINTS, dimx );
//    std::ofstream est_trajectory_file(flightDataPath + "estimated_trajectory.txt", std::ios::out);
//
//    if (!est_trajectory_file.fail()) {
//        for (int i = 0; i < trajectory.size1(); i = i + dimx) {
//            std::vector<double> tmp = trajectory(Slice(i, i + dimx), 0).nonzeros();
//            for (uint j = 0; j < tmp.size(); j++) {
//                est_trajectory_file << tmp[j] << " ";
//            }
//            est_trajectory_file << "\n";
//        }
//    }
//    est_trajectory_file.close();

    /** update parameter file */
    YAML::Node config = YAML::LoadFile(kite_params_file);
    config["geometry"]["imu_pitch_offs"] = new_params_vec[0];

    /** LONGITUDINAL IDENTIFICATION PARAMETERS ---------------------------------------------------------------------- */
    if (flight.maneuver == "identPitch") {
        config["aerodynamic"]["CL0"] = new_params_vec[1];
        config["aerodynamic"]["CLa_total"] = new_params_vec[2];

        config["aerodynamic"]["CD0_total"] = new_params_vec[3];
        config["aerodynamic"]["Cm0"] = new_params_vec[4];
        config["aerodynamic"]["Cma"] = new_params_vec[5];

        config["aerodynamic"]["CLq"] = new_params_vec[6];
        config["aerodynamic"]["Cmq"] = new_params_vec[7];

        config["aerodynamic"]["CLde"] = new_params_vec[8];
        config["aerodynamic"]["Cmde"] = new_params_vec[9];
    }
    /** END OF LONGITUDINAL IDENTIFICATION PARAMETERS --------------------------------------------------------------- */

    /** LATERAL IDENTIFICATION PARAMETERS --------------------------------------------------------------------------- */
    if (flight.maneuver == "identRoll") {
        config["aerodynamic"]["CYb"] = new_params_vec[1];
        config["aerodynamic"]["Cn0"] = new_params_vec[2];
        config["aerodynamic"]["Cnb"] = new_params_vec[3];
        config["aerodynamic"]["Clb"] = new_params_vec[4];

        config["aerodynamic"]["CYr"] = new_params_vec[5];
        config["aerodynamic"]["Cnr"] = new_params_vec[6];
        config["aerodynamic"]["Clr"] = new_params_vec[7];
        config["aerodynamic"]["CYp"] = new_params_vec[8];
        config["aerodynamic"]["Clp"] = new_params_vec[9];
        config["aerodynamic"]["Cnp"] = new_params_vec[10];

        config["aerodynamic"]["CYdr"] = new_params_vec[11];
        config["aerodynamic"]["Cndr"] = new_params_vec[12];
        config["aerodynamic"]["Cldr"] = new_params_vec[13];
        config["aerodynamic"]["Clda"] = new_params_vec[14];
        config["aerodynamic"]["Cnda"] = new_params_vec[15];
    }
    /** END OF LATERAL IDENTIFICATION PARAMETERS -------------------------------------------------------------------- */

    // config["tether"]["length"] = new_params_vec[23];
    // config["tether"]["Ks"] = new_params_vec[24];
    // config["tether"]["Kd"] = new_params_vec[25];
    // config["tether"]["rx"] = new_params_vec[26];
    // config["tether"]["rz"] = new_params_vec[27];

    std::ofstream fout(flightDataPath + "easy_glider_4_new.yaml");
    fout << config;
}
