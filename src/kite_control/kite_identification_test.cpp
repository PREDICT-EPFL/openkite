#include "kiteNMPF.h"
#include "integrator.h"
#include <fstream>
#include <utility>
#include "pseudospectral/chebyshev.hpp"

using namespace casadi;

static DM flipDM(const DM &matrix, int axis) {
    switch (axis) {
        case 1:
            return DM::vertcat({matrix(Slice(-1, 0, -1), Slice()), matrix(0, Slice())});
        case 2:
            return DM::horzcat({matrix(Slice(), Slice(-1, 0, -1)), matrix(Slice(), 0)});
    }
}

void set_relative_parameter_bounds(DM &LBP, DM &UBP, const int &index,
                                   const DM &REF_P, const double &LB_factor, const double &UB_factor) {

    LBP(index) = REF_P(index) + LB_factor * fabs(REF_P(index));
    UBP(index) = REF_P(index) + UB_factor * fabs(REF_P(index));

}

void set_absolute_parameter_bounds(DM &LBP, DM &UBP, const int &index,
                                   const double &absoluteLB, const double &absoluteUB) {

    LBP(index) = absoluteLB;
    UBP(index) = absoluteUB;

}

void printSingleParameterOutput(const int paramNumber, const std::string &paramName,
                                const double &paramFound, const double &param_LB, const double &param_UB,
                                const double &param_sens, const double &sensitivity_max,
                                const double &lagrangeScore = 0) {

    const int p = 4;
    const int w = 1 + 2 + 1 + p; // minus, 2 digits, decimal, 4 digits

    /* Write line */
    std::cout << std::fixed << std::setprecision(p);
    std::cout << std::left;

    std::cout << "(" << std::setw(2) << std::right << paramNumber << std::left << ")"
              << " "
              << std::setw(9) << paramName
              << " = "
              << std::setw(w) << std::right << paramFound << std::left;

    { /* Bound visualization ---------------------------------------------------------------------------------------- */

        /* Lower bound */
        std::cout << std::setw(4) << " "
                  << std::setw(w) << std::right << param_LB
                  << " |";

        { /* Body */

            int nBins = 40;
            double percentagePositionWithinBounds = (paramFound - param_LB) / (param_UB - param_LB);
            int intPositionWithinBounds = static_cast<int>(std::round(percentagePositionWithinBounds * nBins - 0.5));

            if (intPositionWithinBounds <= -1) intPositionWithinBounds = 0;
            if (intPositionWithinBounds >= nBins) intPositionWithinBounds = nBins - 1;
            for (int i = 0; i < nBins; ++i) {
                if (i == intPositionWithinBounds)
                    std::cout << "X";
                else if (i == (nBins / 2))
                    std::cout << "|";
                else
                    std::cout << "-";
            }
        }

        /* Upper bound */
        std::cout << "| " << std::setw(w) << param_UB;
    }

    std::cout << std::setw(4) << "";

    { /* Sensitivity visualization ---------------------------------------------------------------------------------- */

        /* Lower bound */
        std::cout << std::setw(w) << std::right << "0"
                  << " |";

        /* Body */
        int nBins = 40;
        double percentagePositionWithinBounds = param_sens / sensitivity_max;
        int intPositionWithinBounds = static_cast<int>(std::round(percentagePositionWithinBounds * nBins - 0.5));

        if (intPositionWithinBounds <= -1) intPositionWithinBounds = 0;
        if (intPositionWithinBounds >= nBins) intPositionWithinBounds = nBins - 1;
        for (int i = 0; i < nBins; ++i) {
            if (i == intPositionWithinBounds)
                std::cout << "X";
            else
                std::cout << "-";
        }

        /* Upper bound */
        std::cout << "> " << std::setw(w) << param_sens;


    }
    std::cout << "\n";
}

int main(void) {

    /// 5 lines to adapt by the user ///

    struct FlightMetaData {
        int session{0};
        int number{0};
        int seq{0};

        /// 1. Identification maneuver ///
        //std::string maneuver = "identRoll";
        std::string maneuver = "identPitch";

        std::string resampleMethod = "cheb";
    };

    /// 2. ///
    /** COMMENT / UNCOMMENT THE SUITING PARAMETERS IN KITE.CPP ! **/

    /// 3. lon: 10, lat: 14 identification parameters ///
    const int dimp = 10;

    /// 4. Should be constant for sequences of the same maneuver. Get numbers from seqInfo.txt! ///
    const int DATA_POINTS = 346;
    const int poly_order = 3;
    const int num_segments = 115;

    /// 5. ///
    std::string flightDataDir = "/home/johannes/identification/processed_flight_data/";

    /** Load identification schedule **/
    std::vector<FlightMetaData> flights;
    {
        std::string id_schedule_maneuver;
        int id_schedule_session{9};
        int id_schedule_flightNumber{2};
        int id_schedule_seq{2};

        std::ifstream id_schedule_file(flightDataDir + "identPitch_identSchedule.txt", std::ios::in);
        if (!id_schedule_file.fail()) {

            while (id_schedule_file >> id_schedule_session &&
                   id_schedule_file >> id_schedule_flightNumber &&
                   id_schedule_file >> id_schedule_seq) {

                FlightMetaData flightToSchedule;
                flightToSchedule.session = id_schedule_session;
                flightToSchedule.number = id_schedule_flightNumber;
                flightToSchedule.seq = id_schedule_seq;

                flights.push_back(flightToSchedule);
            }

        } else {
            std::cout << "Could not open : id schedule file \n";
            id_schedule_file.clear();
        }
    }

//    struct FlightMetaData {
//        /// 1. ///
//        int session = 9;
//        int number = 2;
//        int seq = 2;
//
//        /// 2. ///
//        //std::string maneuver = "identRoll";
//        std::string maneuver = "identPitch";
//
//        /// 3. ///
//        std::string resampleMethod = "cheb";
//        //std::string resampleMethod = "equal";
//    } flight;

    /** Sequence loop ---------------------------------------------------------------------------------------------- **/
    for (const auto &flight : flights) {

        std::cout << std::left
                  << "\n\n"
                  << "\n" << std::setw(8) << "Session" << " " << flight.session
                  << "\n" << std::setw(8) << "Flight" << " " << flight.number
                  << "\n" << std::setw(8) << "Sequence" << " " << flight.seq
                  << "\n";

        std::string flightDataPath = flightDataDir
                                     + "session_" + std::to_string(flight.session)
                                     + "/flight_" + std::to_string(flight.number);
        if (!flight.maneuver.empty()) flightDataPath.append("/" + flight.maneuver);
        flightDataPath.append("/seq_" + std::to_string(flight.seq) + "/");

        std::cout << "Flight data path: " << flightDataPath << "\n";

        /// Param In file  and number of iterations is set by identConfig file. ///
        std::string kite_params_in_dir;
        std::string kite_params_in_filename;
        int nIt{0};
        std::string kite_params_out_file;

        /* Read in identConfig file */
        std::ifstream id_config_file(flightDataPath + "identConfig.txt", std::ios::in);
        if (!id_config_file.fail()) {
            id_config_file >> kite_params_in_dir;
            id_config_file >> kite_params_in_filename;
            id_config_file >> nIt;

            kite_params_in_dir.append("/");

        } else {
            std::cout << "Could not open : identConfig file \n";
            id_config_file.clear();
        }

        /** Ident iteration loop ----------------------------------------------------------------------------------- **/
        for (int iIt = 0; iIt < nIt; ++iIt) {
            std::cout << "\nIteration " << iIt + 1 << "\n";

            std::string kite_params_in_file;

            if (iIt == 0) {
                kite_params_in_file = kite_params_in_dir + kite_params_in_filename + ".yaml";
            } else {
                /* There was an iteration before, use iterated param file. */
                kite_params_in_file = kite_params_out_file;
            }

            std::cout << "Kite param file IN: " << kite_params_in_file << "\n";

            /** define kite dynamics */
            KiteProperties kite_props = kite_utils::LoadProperties(kite_params_in_file);

            /** Load wind data **/
            std::ifstream id_wind_file(flightDataPath + "seqInfo.txt", std::ios::in);
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
            //std::cout << "Measured state trajectory: size " << id_data.size() << id_data << "\n";

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
            //std::cout << "Measured Control trajectory: size " << id_control.size() << id_control << "\n";


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
            /* Parameter struct for sorting by sensitivity in descending order */
            struct Parameter {
                const int id{};
                const std::string name;
                double sensitivity{};

                Parameter() = default;

                Parameter(const int &id_,
                          std::string name_) :

                        id(id_),
                        name(std::move(name_)) {}

                bool operator<(const Parameter &param) const {
                    return sensitivity > param.sensitivity;
                }
            };

            std::list<Parameter> parameterList;

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
            double imuPitchOffset = kite_props.Geometry.ImuPitchOffset_deg * M_PI / 180.0;

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
                                     imuPitchOffset,
                                     CL0, CLa_tot,
                                     CD0_tot, Cm0, Cma,
                                     CLq, Cmq,
                                     CLde, Cmde,
                                    }); // 10 longitudinal parameters

                /** parameter bounds */
                LBP = REF_P;
                UBP = REF_P;

                set_absolute_parameter_bounds(LBP, UBP, b + 0, -10.0 * M_PI / 180.0,
                                              10.0 * M_PI / 180.0); // imuPitchOffset
                parameterList.emplace_back(0, "PitchOffs");

                set_relative_parameter_bounds(LBP, UBP, b + 1, REF_P, -0.1, 0.1); // CL0
                set_relative_parameter_bounds(LBP, UBP, b + 2, REF_P, -0.05, 0.1); // CLa_tot
                parameterList.emplace_back(1, "CL0");
                parameterList.emplace_back(2, "CLa_tot");

                set_relative_parameter_bounds(LBP, UBP, b + 3, REF_P, -0.1, 0.25); // CD0_tot
                set_relative_parameter_bounds(LBP, UBP, b + 4, REF_P, -0.5, 0.5); // Cm0
                set_relative_parameter_bounds(LBP, UBP, b + 5, REF_P, -0.1, 0.3); // Cma
                parameterList.emplace_back(3, "CD0_tot");
                parameterList.emplace_back(4, "Cm0");
                parameterList.emplace_back(5, "Cma");

                set_relative_parameter_bounds(LBP, UBP, b + 6, REF_P, -0.2, 0.2); // CLq
                set_relative_parameter_bounds(LBP, UBP, b + 7, REF_P, -0.3, 0.3); // Cmq
                parameterList.emplace_back(6, "CLq");
                parameterList.emplace_back(7, "Cmq");

                set_relative_parameter_bounds(LBP, UBP, b + 8, REF_P, -0.5, 0.5); // CLde
                set_relative_parameter_bounds(LBP, UBP, b + 9, REF_P, -0.5, 0.5); // Cmde
                parameterList.emplace_back(8, "CLde");
                parameterList.emplace_back(9, "Cmde");

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
                                     CYb, Cnb, Clb,
                                     CYr, Cnr, Clr, CYp, Clp, Cnp,
                                     CYdr, Cndr, Cldr, Clda, Cnda,
                                    }); // 14 lateral parameters

                /** parameter bounds */
                LBP = REF_P;
                UBP = REF_P;

                set_relative_parameter_bounds(LBP, UBP, b + 0, REF_P, -0.5, 0.5);  // CYb
                set_relative_parameter_bounds(LBP, UBP, b + 1, REF_P, -0.5, 0.5);  // Cnb
                set_relative_parameter_bounds(LBP, UBP, b + 2, REF_P, -0.5, 0.5);  // Clb
                parameterList.emplace_back(0, "CYb");
                parameterList.emplace_back(1, "Cnb");
                parameterList.emplace_back(2, "Clb");

                set_relative_parameter_bounds(LBP, UBP, b + 3, REF_P, -0.3, 0.3);  // CYr
                set_relative_parameter_bounds(LBP, UBP, b + 4, REF_P, -0.5, 0.5);  // Cnr
                set_relative_parameter_bounds(LBP, UBP, b + 5, REF_P, -0.5, 0.5);  // Clr
                set_relative_parameter_bounds(LBP, UBP, b + 6, REF_P, -0.5, 0.5);  // CYp
                set_relative_parameter_bounds(LBP, UBP, b + 7, REF_P, -0.5, 0.5);  // Clp
                set_relative_parameter_bounds(LBP, UBP, b + 8, REF_P, -0.3, 1.0);  // Cnp
                parameterList.emplace_back(3, "CYr");
                parameterList.emplace_back(4, "Cnr");
                parameterList.emplace_back(5, "Clr");
                parameterList.emplace_back(6, "CYp");
                parameterList.emplace_back(7, "Clp");
                parameterList.emplace_back(8, "Cnp");

                set_relative_parameter_bounds(LBP, UBP, b + 9, REF_P, -0.5, 0.5);  // CYdr
                set_relative_parameter_bounds(LBP, UBP, b + 10, REF_P, -0.5, 0.5);  // Cndr
                set_relative_parameter_bounds(LBP, UBP, b + 11, REF_P, -0.5, 0.5);  // Cldr
                set_relative_parameter_bounds(LBP, UBP, b + 12, REF_P, -0.5, 0.5);  // Clda
                set_relative_parameter_bounds(LBP, UBP, b + 13, REF_P, -0.5, 0.5);  // Cnda
                parameterList.emplace_back(9, "CYdr");
                parameterList.emplace_back(10, "Cndr");
                parameterList.emplace_back(11, "Cldr");
                parameterList.emplace_back(12, "Clda");
                parameterList.emplace_back(13, "Cnda");

                Q = SX::diag(SX({0, 1e2, 0,
                                 1e2, 0, 1e2,
                                 1e2, 1e2, 1e1,
                                 1e2, 1e2, 1e2, 1e2}));
            }
            /** END OF LATERAL IDENTIFICATION PARAMETERS -------------------------------------------------------------------- */

            std::cout << "OK: Parameter preparation\n";

            /** ----------------------------------------------------------------------------------*/
            double tf = static_cast<double>(id_data(0, id_data.size2() - 1) - id_data(0, 0));

            std::cout << "tf = " << tf << "\n";

            Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> spectral;
            SX diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, tf);
            diff_constr = diff_constr(casadi::Slice(0, num_segments * poly_order * dimx));

            /* Vectors, length: states * colocation points */
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

            /* Rows: States, Columns: Time steps */
            SX varx_ = SX::reshape(varx, dimx, DATA_POINTS);
            //std::cout << "varx_:" << varx_ << "\n";

            SX fitting_error = 0;
            for (uint j = 0; j < DATA_POINTS; ++j) {
                SX measurement = id_data(Slice(1, id_data.size1()), j);
                SX error = measurement - varx_(Slice(0, varx_.size1()), varx_.size2() - j - 1);
                fitting_error += static_cast<double>(1.0 / DATA_POINTS) * SX::sum1(SX::mtimes(Q, pow(error, 2)));
            }
            Function fitting_error_func = Function("fitting_error", {varx_}, {fitting_error});

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

            /** update parameter file */
            YAML::Node kite_params_yaml = YAML::LoadFile(kite_params_in_file);

            /** LONGITUDINAL IDENTIFICATION PARAMETERS ---------------------------------------------------------------------- */
            if (flight.maneuver == "identPitch") {
                kite_params_yaml["geometry"]["imu_pitch_offs_deg"] = new_params_vec[0] * 180.0 / M_PI;

                kite_params_yaml["aerodynamic"]["CL0"] = new_params_vec[1];
                kite_params_yaml["aerodynamic"]["CLa_total"] = new_params_vec[2];

                kite_params_yaml["aerodynamic"]["CD0_total"] = new_params_vec[3];
                kite_params_yaml["aerodynamic"]["Cm0"] = new_params_vec[4];
                kite_params_yaml["aerodynamic"]["Cma"] = new_params_vec[5];

                kite_params_yaml["aerodynamic"]["CLq"] = new_params_vec[6];
                kite_params_yaml["aerodynamic"]["Cmq"] = new_params_vec[7];

                kite_params_yaml["aerodynamic"]["CLde"] = new_params_vec[8];
                kite_params_yaml["aerodynamic"]["Cmde"] = new_params_vec[9];
            }
            /** END OF LONGITUDINAL IDENTIFICATION PARAMETERS --------------------------------------------------------------- */

            /** LATERAL IDENTIFICATION PARAMETERS --------------------------------------------------------------------------- */
            if (flight.maneuver == "identRoll") {
                kite_params_yaml["aerodynamic"]["CYb"] = new_params_vec[0];
                kite_params_yaml["aerodynamic"]["Cn0"] = new_params_vec[1];
                kite_params_yaml["aerodynamic"]["Cnb"] = new_params_vec[2];
                kite_params_yaml["aerodynamic"]["Clb"] = new_params_vec[3];

                kite_params_yaml["aerodynamic"]["CYr"] = new_params_vec[4];
                kite_params_yaml["aerodynamic"]["Cnr"] = new_params_vec[5];
                kite_params_yaml["aerodynamic"]["Clr"] = new_params_vec[6];
                kite_params_yaml["aerodynamic"]["CYp"] = new_params_vec[7];
                kite_params_yaml["aerodynamic"]["Clp"] = new_params_vec[8];
                kite_params_yaml["aerodynamic"]["Cnp"] = new_params_vec[9];

                kite_params_yaml["aerodynamic"]["CYdr"] = new_params_vec[10];
                kite_params_yaml["aerodynamic"]["Cndr"] = new_params_vec[11];
                kite_params_yaml["aerodynamic"]["Cldr"] = new_params_vec[12];
                kite_params_yaml["aerodynamic"]["Clda"] = new_params_vec[13];
                kite_params_yaml["aerodynamic"]["Cnda"] = new_params_vec[14];
            }
            /** END OF LATERAL IDENTIFICATION PARAMETERS -------------------------------------------------------------------- */

            // config["tether"]["length"] = new_params_vec[23];
            // config["tether"]["Ks"] = new_params_vec[24];
            // config["tether"]["Kd"] = new_params_vec[25];
            // config["tether"]["rx"] = new_params_vec[26];
            // config["tether"]["rz"] = new_params_vec[27];

            /** Save new parameters to file */
            std::string kite_params_out_dir = flightDataPath;
            kite_params_out_dir += "params_new/";
            kite_params_out_dir += kite_params_in_filename + "/";

            kite_params_out_file = kite_params_out_dir + "afterIt_" + std::to_string(iIt + 1) + ".yaml";
            std::ofstream fout(kite_params_out_file);
            fout << kite_params_yaml;


            /** Save estimated trajectory */
            std::ofstream est_trajectory_file(kite_params_out_dir + "afterIt_" + std::to_string(iIt + 1) +
                                              "_estimatedTrajectory.txt", std::ios::out);

            DM trajectory_vect = result(Slice(0, varx.size1()));

            /* Rows: States, Columns: Time steps */
            DM estimated_trajectory_woTime = flipDM(DM::reshape(trajectory_vect, dimx, DATA_POINTS), 2);

            /* Add original times as first row */
            DM estimated_trajectory = DM::vertcat({id_data(0, Slice(0, id_data.size2())), estimated_trajectory_woTime});

            std::ofstream trajectory_file(kite_params_out_dir + "afterIt_" + std::to_string(iIt + 1) +
                                          "_estimatedTrajectory.txt", std::ios::out);
            trajectory_file.precision(15);
            for (int iTimeStep = 0; iTimeStep < estimated_trajectory.size2(); ++iTimeStep) {
                for (int jState = 0; jState < estimated_trajectory.size1(); ++jState) {
                    trajectory_file << static_cast<double>(estimated_trajectory(jState, iTimeStep));
                    if (jState < estimated_trajectory.size1() - 1) { trajectory_file << " "; }
                }
                trajectory_file << "\n";
            }
            trajectory_file.close();

            /** Calculate fitting error (cost function) */
            DM fitting_error_evaluated = fitting_error_func(DMVector{estimated_trajectory_woTime})[0];

            /** Visualize parameters found within their bounds */
            /* Add sensibility values to parameter list (order of coding) */
            for (Parameter &param : parameterList) {
                param.sensitivity = std::abs(static_cast<double>(param_sens(param.id)));
            }

            /* Sort parameter list by sensitivity */
            parameterList.sort();

            double sensitivity_max = parameterList.begin()->sensitivity;

            /* Header */
            std::cout << "\n\n";
            std::cout << std::left;
            std::cout << std::setw(4) << "No."
                      << " "
                      << std::setw(9) << "Name"
                      << "   "
                      << std::setw(8) << "Value"

                      << std::setw(4) << ""

                      << std::setw(8) << "Lw.Bound"
                      << "  "
                      << std::setw(40) << "Value position in bounds"
                      << "  "
                      << std::setw(8) << "Up.Bound"

                      << std::setw(4) << ""

                      << std::setw(8) << ""
                      << "  "
                      << std::setw(40) << "Sensitivity"
                      << "\n";

            /* Print parameter outputs */
            for (Parameter &param : parameterList) {
                int i = param.id;

                printSingleParameterOutput(i, param.name,
                                           new_params_vec[i], static_cast<double>(LBP(i)), static_cast<double>(UBP(i)),
                                           param.sensitivity, sensitivity_max);
            }
            std::cout << "\nFitting error = " << static_cast<double>(fitting_error_evaluated) << "\n";

            /** Save paramInfo file **/
            std::ofstream param_diag_file(kite_params_out_dir + "afterIt_" + std::to_string(iIt + 1) + "_paramInfo.txt",
                                          std::ios::out);

            param_diag_file << "ParamId" << " "
                            << "ParamName" << " "
                            << "ParamValue" << " "
                            << "LwBound" << " "
                            << "UpBound" << " "
                            << "Sensitivity" << " "
                            << "FittingError" << "\n";

            for (Parameter &param : parameterList) {
                int i = param.id;

                param_diag_file << std::fixed << std::setprecision(15);

                param_diag_file << i << " "
                                << param.name << " "
                                << new_params_vec[i] << " "
                                << static_cast<double>(LBP(i)) << " "
                                << static_cast<double>(UBP(i)) << " "
                                << param.sensitivity << " "
                                << static_cast<double>(fitting_error_evaluated) << "\n";

            }
            param_diag_file.close();

        } /** End of ident iteration loop -------------------------------------------------------------------------- **/

    } /** End of sequence loop ------------------------------------------------------------------------------------- **/

}
