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


struct OptProblemProperties {

    const int dimx;
    const int dimu;
    const int dimp;

    const int DATA_POINTS;
    const int poly_order;
    const int num_segments;

    OptProblemProperties(const int &dimx_, const int &dimu_, const int &dimp_,
                         const int &DATA_POINTS_, const int &poly_order_, const int &num_segments_) :

            dimx(dimx_), dimu(dimu_), dimp(dimp_),
            DATA_POINTS(DATA_POINTS_), poly_order(poly_order_), num_segments(num_segments_) {}

};

/* Data in */
struct FlightMetaData {
    int session{0};
    int number{0};
    int seq{0};

    std::string resampleMethod = "cheb";
};

/* Parameter struct for sorting by sensitivity in descending order */
struct Parameter {

    enum BoundType {
        REL,
        ABS
    };

    const int id{};
    const std::string name;
    const double refValue{};
    const double lowerBound{};
    const double upperBound{};

    double sensitivity{};

    Parameter() = default;

    /* Initialization with REL/ABS bound option */
    Parameter(const int &id_,
              std::string &name_,
              const double &refValue_,
              const BoundType boundType_,
              const double &lowerBound_,
              const double &upperBound_) :

            id(id_),
            name(std::move(name_)),
            refValue(refValue_),
            lowerBound((boundType_ == REL) ? (refValue_ - lowerBound_ * std::abs(refValue_)) : lowerBound_),
            upperBound((boundType_ == REL) ? (refValue_ + upperBound_ * std::abs(refValue_)) : upperBound_) {}

    /* Initialization with REL/ABS bound option */
    Parameter(const int &id_,
              std::string &name_,
              const double &refValue_,
              const double &absoluteLowerBound_,
              const double &absoluteUpperBound_) :

            id(id_),
            name(std::move(name_)),
            refValue(refValue_),
            lowerBound(absoluteLowerBound_),
            upperBound(absoluteUpperBound_) {}

    bool operator<(const Parameter &param) const {
        return sensitivity > param.sensitivity;
    }
};

std::vector<FlightMetaData> readIn_identSchedule(const std::string &filedir, const std::string &identManeuver) {


    std::string filepath = filedir + identManeuver + "_identSchedule.txt";

    std::vector<FlightMetaData> flights;

    int id_schedule_session{9};
    int id_schedule_flightNumber{2};
    int id_schedule_seq{2};

    std::ifstream id_schedule_file(filepath, std::ios::in);
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

    return flights;
}

void get_seqDir(const std::string &flightDataDir, const std::string &identManeuver, const FlightMetaData &flight,

                std::string &seq_dir) {

    seq_dir = flightDataDir
              + "session_" + std::to_string(flight.session)
              + "/flight_" + std::to_string(flight.number);
    if (!identManeuver.empty()) seq_dir.append("/" + identManeuver);
    seq_dir.append("/seq_" + std::to_string(flight.seq) + "/");
}


void readIn_identificationData(const std::string &filepath, const int &dimension, const int &DATA_POINTS,

                               DM &identData) {

    /* File to read in:
     * Rows: Timesteps, Columns: States/Controls
     *
     * Each row contains a vector of the states at that time. */

    /** Load identification data */
    std::ifstream id_data_file(filepath, std::ios::in);

    /* Columns: Timesteps, Rows: States/Controls
     * Each column contains a vector of the states at that time */
    DM id_data = DM::zeros(dimension, DATA_POINTS);

    /** load data */
    if (!id_data_file.fail()) {

        for (uint iDATA_POINT = 0; iDATA_POINT < DATA_POINTS; ++iDATA_POINT) {
            /* Loop through rows */

            for (uint jDim = 0; jDim < dimension; ++jDim) {
                /* Loop through columns */

                double entry;
                id_data_file >> entry;
                id_data(jDim, iDATA_POINT) = entry;
//                identData(jDim, iDATA_POINT) = entry;

            }
        }
    } else {
        std::cout << "Could not open id data file.\n";
        id_data_file.clear();
    }

    identData = id_data;
}

void readIn_sequenceInfo(const std::string &filepath,

                         int &windFrom_deg, double &windSpeed, double &tf) {

    /** Load wind data **/
    std::ifstream id_seqInfo_file(filepath, std::ios::in);
    if (!id_seqInfo_file.fail()) {

        std::string keyBuffer;

        id_seqInfo_file >> windFrom_deg;
        id_seqInfo_file >> windSpeed;

        id_seqInfo_file >> keyBuffer;
        id_seqInfo_file >> tf;

    } else {
        std::cout << "Could not open: Sequence info file\n";
        id_seqInfo_file.clear();
    }
}
//double readIn_sequenceInfo(const std::string &filepath) {
//
//    int windFrom_deg{0};
//    double windSpeed{0};
//    double tf{0};
//    readIn_sequenceInfo(filepath, windFrom_deg, windSpeed, tf);
//
//    return tf;
//}
//void readIn_sequenceInfo(const std::string &filepath,
//
//                         int &windFrom_deg, double &windSpeed) {
//
//    double tf{0};
//    readIn_sequenceInfo(filepath,
//
//                        windFrom_deg, windSpeed, tf);
//}

void readIn_identConfig(const std::string &filepath,

                        std::string &kite_params_in_dir, std::string &kite_params_in_filename,
                        int &nIt, bool &kite_params_warmStart) {

    int kite_params_warmStart_int{0};

    /* Read in identConfig file */
    std::ifstream id_config_file(filepath, std::ios::in);
    if (!id_config_file.fail()) {
        id_config_file >> kite_params_in_dir;
        id_config_file >> kite_params_in_filename;
        id_config_file >> nIt;
        id_config_file >> kite_params_warmStart; //_int;

        //kite_params_warmStart = static_cast<bool>(kite_params_warmStart_int);

        kite_params_in_dir.append("/");

    } else {
        std::cout << "Could not open : identConfig file \n";
        id_config_file.clear();
    }
}


void readIn_sequenceData(const std::string &seq_dir, const FlightMetaData &flight,
                         const int &dimx, const int &dimu, const int &DATA_POINTS,

                         DM &id_time, DM &id_states_woTime, DM &id_controls_rev_woTime,
                         int &windFrom_deg, double &windSpeed, double &tf,
                         std::string &kite_baseParams_dir, std::string &kite_baseParams_filename,
                         int &nIt, bool &kite_params_warmStart) {

    /* States */
    DM id_states_wTime;
    readIn_identificationData(seq_dir + flight.resampleMethod + "_states.txt", 1 + dimx, DATA_POINTS,

                              id_states_wTime);

    id_time = id_states_wTime(0, Slice(0, DATA_POINTS));
    id_states_woTime = id_states_wTime(Slice(1, 1 + dimx), Slice(0, DATA_POINTS));

    /* Controls */
    DM id_controls_wTime;
    readIn_identificationData(seq_dir + flight.resampleMethod + "_controls.txt", 1 + dimu, DATA_POINTS,

                              id_controls_wTime);

    DM id_controls_rev = flipDM(id_controls_wTime, 2);
    id_controls_rev_woTime = id_controls_rev(Slice(1, 1 + dimu), Slice(0, DATA_POINTS));

    /* Sequence info */
    readIn_sequenceInfo(seq_dir + "seqInfo.txt",

                        windFrom_deg, windSpeed, tf);

    /* Ident config */
    readIn_identConfig(seq_dir + "identConfig.txt",

                       kite_baseParams_dir, kite_baseParams_filename, nIt, kite_params_warmStart);
}


/* Set up */
void get_costMatrix(const KiteDynamics::IdentMode identMode,

                    DM &Q) {

    if (identMode == KiteDynamics::IdentMode::LONGITUDINAL) {

        Q = SX::diag(SX({1e2, 0, 1e2,
                         0, 1e2, 0,
                         0, 0, 1e2,
                         1e2, 0, 1e2, 0}));
    }

    if (identMode == KiteDynamics::IdentMode::LATERAL) {

        Q = SX::diag(SX({1e1, 1e2, 1e1,
                         1e2, 1e1, 1e2,
                         1e2, 1e2, 1e2,
                         1e2, 1e2, 1e2, 1e2}));
    }

}

void get_kiteDynamics(const std::string &filepath, const double windFrom_deg, const double windSpeed,
                      const KiteDynamics::IdentMode &identMode,

                      KiteDynamics &kite, KiteDynamics &kite_int) {

    KiteProperties kite_props = kite_utils::LoadProperties(filepath);
    kite_props.Wind.WindFrom_deg = windFrom_deg;
    kite_props.Wind.WindSpeed = windSpeed;

    AlgorithmProperties algo_props;
    algo_props.Integrator = CVODES;
    algo_props.sampling_time = 0.02;
    //algo_props.sampling_time = static_cast<double>(id_data(0, 1) - id_data(0, 0));

    kite = KiteDynamics(kite_props, algo_props, identMode);
    kite_int = KiteDynamics(kite_props, algo_props); //integration model
}

void create_nlpSolver(const SX &opt_var, const SX &fitting_error, const SX &diff_constr,

                      Function &nlpSolver) {

    /** formulate NLP */
    SXDict NLP;
    NLP["x"] = opt_var;
    NLP["f"] = fitting_error;
    NLP["g"] = diff_constr;

    Dict OPTS;
    OPTS["ipopt.linear_solver"] = "mumps";
    OPTS["ipopt.print_level"] = 5;
    OPTS["ipopt.tol"] = 1e-4;
    OPTS["ipopt.acceptable_tol"] = 1e-4;
    OPTS["ipopt.warm_start_init_point"] = "yes";
    //OPTS["ipopt.max_iter"]       = 20;

    nlpSolver = nlpsol("solver", "ipopt", NLP, OPTS);
}

void setup_optimizationParameters(const KiteProperties &kiteProps, const std::string &identManeuver,

                                  DM &REF_P, DM &LBP, DM &UBP, std::list<Parameter> &parameterList) {

    /** enviromental constants */
    double g = 9.80665; /** gravitational acceleration [m/s2] [WGS84] */
    double rho = 1.2985; /** standard atmospheric density [kg/m3] [standard Atmosphere 1976] */

    /** --------------------- **/
    /** Wind properties       **/
    /** --------------------- **/
    double windFrom_deg = kiteProps.Wind.WindFrom_deg;
    double windSpeed = kiteProps.Wind.WindSpeed;

    /** --------------------- **/
    /** Geometric parameters  **/
    /** --------------------- **/
    double imuPitchOffset_deg = kiteProps.Geometry.ImuPitchOffset_deg;

    double b = kiteProps.Geometry.WingSpan;
    double c = kiteProps.Geometry.MAC;
    double AR = kiteProps.Geometry.AspectRatio;
    double S = kiteProps.Geometry.WingSurfaceArea;
    //double lam = KiteProps.Geometry.TaperRatio;
    //double St = KiteProps.Geometry.HTailsurface;
    //double lt = KiteProps.Geometry.TailLeverArm;
    //double Sf = KiteProps.Geometry.FinSurfaceArea;
    //double lf = KiteProps.Geometry.FinLeverArm;
    //double Xac = KiteProps.Geometry.AerodynamicCenter;
    /** @todo: get rid of all magic numbers **/
    //double Xcg = 0.031/c;               /** Center of Gravity (CoG) wrt leading edge [1/c] **/
    //double Vf = (Sf * lf) / (S * b);    /** fin volume coefficient []                      **/
    //double Vh = (lt * St) / (S * c);    /** horizontal tail volume coefficient []          **/

    /** --------------------------- **/
    /** Mass and inertia parameters **/
    /** --------------------------- **/
    double Mass = kiteProps.Inertia.Mass;
    double Ixx = kiteProps.Inertia.Ixx;
    double Iyy = kiteProps.Inertia.Iyy;
    double Izz = kiteProps.Inertia.Izz;
    double Ixz = kiteProps.Inertia.Ixz;

    /** ------------------------------- **/
    /** Aerodynamic parameters          **/
    /** ------------------------------- **/
    double e_o = kiteProps.Aerodynamics.e_oswald;
    double CD0 = kiteProps.Aerodynamics.CD0;


    /* AOA */
    double CL0 = kiteProps.Aerodynamics.CL0;
    double CLa = kiteProps.Aerodynamics.CLa;

    double Cm0 = kiteProps.Aerodynamics.Cm0;
    double Cma = kiteProps.Aerodynamics.Cma;

    /* Sideslip */
    double CYb = kiteProps.Aerodynamics.CYb;

    double Cl0 = kiteProps.Aerodynamics.Cl0;
    double Clb = kiteProps.Aerodynamics.Clb;

    double Cn0 = kiteProps.Aerodynamics.Cn0;
    double Cnb = kiteProps.Aerodynamics.Cnb;

    /* Pitchrate */
    double CLq = kiteProps.Aerodynamics.CLq;
    double Cmq = kiteProps.Aerodynamics.Cmq;

    /* Rollrate */
    double CYp = kiteProps.Aerodynamics.CYp;
    double Clp = kiteProps.Aerodynamics.Clp;
    double Cnp = kiteProps.Aerodynamics.Cnp;

    /* Yawrate */
    double CYr = kiteProps.Aerodynamics.CYr;
    double Clr = kiteProps.Aerodynamics.Clr;
    double Cnr = kiteProps.Aerodynamics.Cnr;


    /** ------------------------------ **/
    /** Aerodynamic effects of control **/
    /** ------------------------------ **/
    /* Elevator */
    double CLde = kiteProps.Aerodynamics.CLde;
    double Cmde = kiteProps.Aerodynamics.Cmde;

    /* Ailerons */
    double Clda = kiteProps.Aerodynamics.Clda;
    double Cnda = kiteProps.Aerodynamics.Cnda;

    /* Rudder */
    double CYdr = kiteProps.Aerodynamics.CYdr;
    double Cldr = kiteProps.Aerodynamics.Cldr;
    double Cndr = kiteProps.Aerodynamics.Cndr;

    //double CL_daoa = -2 * CLa_t * Vh * dw;
    //double Cm_daoa = -2 * CLa_t * Vh * (lt/c) * dw;

    /** ------------------------------ **/
    /**        Tether parameters       **/
    /** ------------------------------ **/
    //    double Lt = KiteProps.Tether.length;
    //    double Ks = KiteProps.Tether.Ks;
    //    double Kd = KiteProps.Tether.Kd;
    //    double rx = KiteProps.Tether.rx;
    //    double ry = KiteProps.Tether.ry;
    //    double rz = KiteProps.Tether.rz;

    //DM REF_P = DM::vertcat({Lt, Ks, Kd, rx, ry, rz});

    /** LONGITUDINAL IDENTIFICATION PARAMETERS ---------------------------------------------------------------------- */
    if (identManeuver == "identPitch") {
        int b = REF_P.size1(); // parameter size before

        REF_P = DM::vertcat({REF_P,
                             imuPitchOffset_deg,

                             CD0,

                             CL0,
                             CLa,

                             Cm0,
                             Cma,

                             CLq,
                             Cmq,

                             CLde,
                             Cmde
                            }); // 10 longitudinal parameters

        /** parameter bounds */
        LBP = UBP = REF_P;

        //set_abs_bounds(b + 0, kiteProps.Geometry.ImuPitchOffset_deg, -0.1, 0.1, REF_P, LBP, UBP);

        set_absolute_parameter_bounds(LBP, UBP, b + 0, -0.1, 0.1); // imuPitchOffset_deg
        parameterList.emplace_back(0, "PitchOffs");

        set_relative_parameter_bounds(LBP, UBP, b + 1, REF_P, -0.5, 0.7); // CD0
        parameterList.emplace_back(1, "CD0");

        set_relative_parameter_bounds(LBP, UBP, b + 2, REF_P, -0.1, 0.5); // CL0
        set_relative_parameter_bounds(LBP, UBP, b + 3, REF_P, -0.05, 0.1); // CLa_
        parameterList.emplace_back(2, "CL0");
        parameterList.emplace_back(3, "CLa");

        set_relative_parameter_bounds(LBP, UBP, b + 4, REF_P, -0.5, 0.5); // Cm0
        set_relative_parameter_bounds(LBP, UBP, b + 5, REF_P, -0.1, 0.3); // Cma
        parameterList.emplace_back(4, "Cm0");
        parameterList.emplace_back(5, "Cma");

        set_relative_parameter_bounds(LBP, UBP, b + 6, REF_P, -0.5, 0.5); // CLq
        set_relative_parameter_bounds(LBP, UBP, b + 7, REF_P, -0.3, 0.3); // Cmq
        parameterList.emplace_back(6, "CLq");
        parameterList.emplace_back(7, "Cmq");

        set_relative_parameter_bounds(LBP, UBP, b + 8, REF_P, -0.5, 0.7); // CLde
        set_relative_parameter_bounds(LBP, UBP, b + 9, REF_P, -0.5, 0.5); // Cmde
        parameterList.emplace_back(8, "CLde");
        parameterList.emplace_back(9, "Cmde");

    }
    /** END OF LONGITUDINAL IDENTIFICATION PARAMETERS --------------------------------------------------------------- */

    /** LATERAL IDENTIFICATION PARAMETERS --------------------------------------------------------------------------- */
    if (identManeuver == "identRoll") {
        int b = REF_P.size1(); // parameter size before

        REF_P = DM::vertcat({REF_P,
                             CYb,

                             Clb,

                             Cnb,

                             CYp,
                             Clp,
                             Cnp,

                             CYr,
                             Clr,
                             Cnr,

                             Clda,
                             Cnda
                            }); // 11 lateral parameters (aileron only)

        /** parameter bounds */
        LBP = UBP = REF_P;

        set_relative_parameter_bounds(LBP, UBP, b + 0, REF_P, -0.5, 0.5);  // CYb
        parameterList.emplace_back(0, "CYb");

        set_relative_parameter_bounds(LBP, UBP, b + 1, REF_P, -0.5, 0.5);  // Clb
        parameterList.emplace_back(1, "Clb");

        set_relative_parameter_bounds(LBP, UBP, b + 2, REF_P, -0.5, 0.5);  // Cnb
        parameterList.emplace_back(2, "Cnb");

        set_relative_parameter_bounds(LBP, UBP, b + 3, REF_P, -0.5, 0.5);  // CYp
        set_relative_parameter_bounds(LBP, UBP, b + 4, REF_P, -0.5, 0.5);  // Clp
        set_relative_parameter_bounds(LBP, UBP, b + 5, REF_P, -0.3, 1.0);  // Cnp
        parameterList.emplace_back(3, "CYp");
        parameterList.emplace_back(4, "Clp");
        parameterList.emplace_back(5, "Cnp");

        set_relative_parameter_bounds(LBP, UBP, b + 6, REF_P, -0.3, 0.3);  // CYr
        set_relative_parameter_bounds(LBP, UBP, b + 7, REF_P, -0.5, 0.5);  // Clr
        set_relative_parameter_bounds(LBP, UBP, b + 8, REF_P, -0.5, 0.5);  // Cnr
        parameterList.emplace_back(6, "CYr");
        parameterList.emplace_back(7, "Clr");
        parameterList.emplace_back(8, "Cnr");

        set_relative_parameter_bounds(LBP, UBP, b + 9, REF_P, -0.5, 0.5);  // Clda
        set_relative_parameter_bounds(LBP, UBP, b + 10, REF_P, -0.5, 0.5);  // Cnda
        parameterList.emplace_back(9, "Clda");
        parameterList.emplace_back(10, "Cnda");

//        set_relative_parameter_bounds(LBP, UBP, b + 11, REF_P, -0.5, 0.5);  // CYdr
//        set_relative_parameter_bounds(LBP, UBP, b + 12, REF_P, -0.5, 0.5);  // Cldr
//        set_relative_parameter_bounds(LBP, UBP, b + 13, REF_P, -0.5, 0.5);  // Cndr
//        parameterList.emplace_back(11, "CYdr");
//        parameterList.emplace_back(12, "Cldr");
//        parameterList.emplace_back(13, "Cndr");

    }
    /** END OF LATERAL IDENTIFICATION PARAMETERS -------------------------------------------------------------------- */

}

void get_optimizationVarBounds(const DM &LBX, const DM &UBX, const int &num_segments, const int &poly_order,
                               const DM &LBU, const DM &UBU,
                               const DM &LBP, const DM &UBP,

                               SX &lbx, SX &ubx) {

    /** Optimization variable bounds (inequality (box) constraints) **/
    /* Add state bounds for each DATA_POINT (replicate vector for one DATA_POINT) */
    lbx = SX::repmat(LBX, num_segments * poly_order + 1, 1);
    ubx = SX::repmat(UBX, num_segments * poly_order + 1, 1);

    /* Add control bounds for each DATA_POINT (entirely from DATA)) */
    lbx = SX::vertcat({lbx, LBU});
    ubx = SX::vertcat({ubx, UBU});

    /* Add parameter bounds (once) */
    lbx = SX::vertcat({lbx, LBP});
    ubx = SX::vertcat({ubx, UBP});
}


/* Data out*/
void write_parameterFile(const std::string &kite_params_in_filepath, const std::string &kite_params_out_filepath,
                         const std::string &identManeuver, const std::vector<double> &new_params_vec) {

    /** update parameter file */
    YAML::Node kite_params_yaml = YAML::LoadFile(kite_params_in_filepath);

    /** LONGITUDINAL IDENTIFICATION PARAMETERS ---------------------------------------------------------------------- */
    if (identManeuver == "identPitch") {

        kite_params_yaml["geometry"]["imu_pitch_offs_deg"] = new_params_vec[0];


        kite_params_yaml["aerodynamic"]["CD0"] = new_params_vec[1];

        /* AOA */
        kite_params_yaml["aerodynamic"]["CL0"] = new_params_vec[2];
        kite_params_yaml["aerodynamic"]["CLa"] = new_params_vec[3];

        kite_params_yaml["aerodynamic"]["Cm0"] = new_params_vec[4];
        kite_params_yaml["aerodynamic"]["Cma"] = new_params_vec[5];


        /* Pitchrate */
        kite_params_yaml["aerodynamic"]["CLq"] = new_params_vec[6];
        kite_params_yaml["aerodynamic"]["Cmq"] = new_params_vec[7];


        /* Aerodynamic effects of control */
        /* Elevator */
        kite_params_yaml["aerodynamic"]["CLde"] = new_params_vec[8];
        kite_params_yaml["aerodynamic"]["Cmde"] = new_params_vec[9];
    }
    /** END OF LONGITUDINAL IDENTIFICATION PARAMETERS --------------------------------------------------------------- */

    /** LATERAL IDENTIFICATION PARAMETERS --------------------------------------------------------------------------- */
    if (identManeuver == "identRoll") {

        /* Sideslip */
        kite_params_yaml["aerodynamic"]["CYb"] = new_params_vec[0];

        kite_params_yaml["aerodynamic"]["Clb"] = new_params_vec[1];

        kite_params_yaml["aerodynamic"]["Cnb"] = new_params_vec[2];


        /* Rollrate */
        kite_params_yaml["aerodynamic"]["CYp"] = new_params_vec[3];
        kite_params_yaml["aerodynamic"]["Clp"] = new_params_vec[4];
        kite_params_yaml["aerodynamic"]["Cnp"] = new_params_vec[5];

        /* Yawrate */
        kite_params_yaml["aerodynamic"]["CYr"] = new_params_vec[6];
        kite_params_yaml["aerodynamic"]["Clr"] = new_params_vec[7];
        kite_params_yaml["aerodynamic"]["Cnr"] = new_params_vec[8];


        /* Aerodynamic effects of control */
        /* Ailerons */
        kite_params_yaml["aerodynamic"]["Clda"] = new_params_vec[9];
        kite_params_yaml["aerodynamic"]["Cnda"] = new_params_vec[10];

//        /* Rudder */
//        kite_params_yaml["aerodynamic"]["CYdr"] = new_params_vec[11];
//        kite_params_yaml["aerodynamic"]["Cldr"] = new_params_vec[12];
//        kite_params_yaml["aerodynamic"]["Cndr"] = new_params_vec[13];
    }
    /** END OF LATERAL IDENTIFICATION PARAMETERS -------------------------------------------------------------------- */

    // config["tether"]["length"] = new_params_vec[23];
    // config["tether"]["Ks"] = new_params_vec[24];
    // config["tether"]["Kd"] = new_params_vec[25];
    // config["tether"]["rx"] = new_params_vec[26];
    // config["tether"]["rz"] = new_params_vec[27];

    std::ofstream fout(kite_params_out_filepath);
    fout << kite_params_yaml;
}

void write_trajectoryFile(const DM &traj, const std::string &filepath) {

    std::ofstream traj_file(filepath, std::ios::out);
    traj_file.precision(15);

    for (int iTimeStep = 0; iTimeStep < traj.size2(); ++iTimeStep) {
        for (int jState = 0; jState < traj.size1(); ++jState) {
            traj_file << static_cast<double>(traj(jState, iTimeStep));
            if (jState < traj.size1() - 1) { traj_file << " "; }
        }
        traj_file << "\n";
    }

    traj_file.close();
}

void print_singleParameterOutput(const int paramNumber, const std::string &paramName,
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

void printWrite_parametersFound(std::list<Parameter> &parameterList, const std::vector<double> &new_params_vec,
                                const DM &LBP, const DM &UBP, const DM &param_sens, const DM &fitting_error_evaluated,
                                const std::string &filepath) {

    //std::cout << "PARAMETER SENSITIVITIES: " << param_sens << "\n";

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

        print_singleParameterOutput(i, param.name,
                                    new_params_vec[i],
                                    static_cast<double>(LBP(i)),
                                    static_cast<double>(UBP(i)),
                                    param.sensitivity, sensitivity_max);
    }
    std::cout << "\nFitting error = " << static_cast<double>(fitting_error_evaluated) << "\n";


    /** Save paramInfo file **/

    std::ofstream param_diag_file(filepath, std::ios::out);

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

}


int main() {

    /// 5 lines to adapt by the user =============================================================================== ///

    std::string flightDataDir = "/home/johannes/identification/processed_flight_data/";

    /// State and Control dimensions ///
    const int dimx = 13;  // v(3) w(3) r(3) q(4)
    const int dimu = 4;   // T elev rud ail

    /// 1. Identification mode ///
    const KiteDynamics::IdentMode identMode = KiteDynamics::IdentMode::LATERAL;

    /// 2. lon: 10, lat: 11 identification parameters ///
    const int dimp = 11;

    /// 3. Should be constant for sequences of the same maneuver. Get numbers from seqInfo.txt! ///
    // pitch / longitudinal
//    const int DATA_POINTS = 346;
//    const int poly_order = 3;
//    const int num_segments = 115;

    // roll / lateral
    const int DATA_POINTS = 325;;
    const int poly_order = 3;
    const int num_segments = 108;

    //OptProblemProperties opp(13, 4, 10, 346, 3, 115);
    //OptProblemProperties opp(dimx, dimu, dimp, DATA_POINTS, poly_order, num_segments);

    const int n_optimizationFailed_max = 3;

    /// ============================================================================================================ ///
    /* User config check */
    if (!(
            (identMode == KiteDynamics::IdentMode::LONGITUDINAL && dimp == 10 && DATA_POINTS == 346)
            ||
            (identMode == KiteDynamics::IdentMode::LATERAL && dimp == 11 && DATA_POINTS == 325))) {

        std::cout << "Identification mode/dimp/DATA_POINTS mismatch. Check user config!\n";
        return EXIT_FAILURE;
    }

    std::string identManeuver;
    if (identMode == KiteDynamics::IdentMode::LONGITUDINAL) {
        identManeuver = "identPitch";
    } else if (identMode == KiteDynamics::IdentMode::LATERAL) {
        identManeuver = "identRoll";
    }

    std::cout << "Running " << identManeuver << " identification of " << dimp << " parameters with "
              << DATA_POINTS << " data points.\n";

    /* Load identification schedule */
    std::vector<FlightMetaData> flights = readIn_identSchedule(flightDataDir, identManeuver);

    /* Cost matrix */
    DM Q;
    get_costMatrix(identMode,

                   Q);

    /* State bounds (constant over all sequences and iterations) */
    DM UBX = DM::vertcat({40, 15, 15,
                          4 * M_PI, 4 * M_PI, 4 * M_PI,
                          300, 300, 0,
                          1.05, 1.05, 1.05, 1.05});

    DM LBX = DM::vertcat({2.0, -15, -15,
                          -4 * M_PI, -4 * M_PI, -4 * M_PI,
                          -300, -300, -300,
                          -1.05, -1.05, -1.05, -1.05});     // for infinity use -DM::inf(1) and DM::inf(1)

    /** Collocation **/
    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> spectral;

    /* Vectors, length: states * colocation points */
    SX varx = spectral.VarX();
    SX varu = spectral.VarU();
    SX varp = spectral.VarP();

    /* Optimization variable contains
     * states at each DATA_POINT
     * controls at each DATA_POINT
     * parameters */
    SX opt_var = SX::vertcat(SXVector{varx, varu, varp});

    /** Sequence loop ============================================================================================== **/
    for (const auto &flight : flights) {

        /* Print current sequence */
        std::cout << std::left
                  << "\n---------------------------------------------------------"
                  << "\nData: "
                  << flight.session << "-" << flight.number << "-" << flight.seq
                  << " (Session-Flight-Sequence)\n";


        /* Construct sequence directory path */
        std::string seq_dir;
        get_seqDir(flightDataDir, identManeuver, flight,

                   seq_dir);
        std::cout << "Directory: " << seq_dir << "\n";


        /** Read in all sequence data **/
        DM id_time;
        DM id_states_woTime;
        DM id_controls_rev_woTime;

        int windFrom_deg{0};
        double windSpeed{0};
        double tf{0};

        std::string kite_baseParams_dir;
        std::string kite_baseParams_filename;
        int nIt{0};
        bool kite_params_warmStart{false};

        readIn_sequenceData(seq_dir, flight, dimx, dimu, DATA_POINTS,

                            id_time, id_states_woTime, id_controls_rev_woTime,
                            windFrom_deg, windSpeed, tf,
                            kite_baseParams_dir, kite_baseParams_filename, nIt, kite_params_warmStart);
        //std::cout << "Measured Control trajectory: size " << id_control.size() << id_control << "\n";
        //std::cout << "tf = " << tf << "\n";


        /** Initial state: First column of id_states **/
        DM id_state0 = id_states_woTime(Slice(0, dimx), 0);

        /** Control bounds **/
        /* (for all DATA_POINTs)
         * lower bound = control value = upper bound, using controls as optimization input */
        DM LBU, UBU;
        LBU = UBU = DM::vec(id_controls_rev_woTime);


        /** Fitting error **/
        SX varx_ = SX::reshape(varx, dimx, DATA_POINTS);
        SX fitting_error = 0;
        for (uint jTimeStep = 0; jTimeStep < DATA_POINTS; ++jTimeStep) {

            /* Measured and to be optimized state at one DATA_POINT */
            SX measured_state = id_states_woTime(Slice(0, dimx), jTimeStep);
            SX state_to_optimize = varx_(Slice(0, dimx), DATA_POINTS - jTimeStep - 1); // time-reverse

            /* Error and sum up */
            SX error = measured_state - state_to_optimize;
            fitting_error += static_cast<double>(1.0 / DATA_POINTS) * SX::sum1(SX::mtimes(Q, pow(error, 2)));

        }
        Function fitting_error_func = Function("fitting_error", {varx_}, {fitting_error});


        /** Kite Dynamics based on baseParams and current sequence's wind conditions **/
        /* Get KiteDynamics based on constant (during optimization) parameters.
         * Parameters to be optimized parameters are read in for each iteration. */
        std::string kite_params_in_filepath = kite_baseParams_dir + kite_baseParams_filename + ".yaml";

        KiteDynamics kite, kite_int;
        get_kiteDynamics(kite_params_in_filepath, windFrom_deg, windSpeed, identMode,

                         kite, kite_int);

        Function DynamicsFunc = kite.getNumericDynamics();
        //SX X = kite.getSymbolicState();
        //SX U = kite.getSymbolicControl();
        SX P = kite.getSymbolicParameters();

        Function ode = kite_int.getNumericDynamics();

        SX diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, tf);
        diff_constr = diff_constr(Slice(0, num_segments * poly_order * dimx));


        /** Create solver **/
        Function nlp_solver_func;
        create_nlpSolver(opt_var, fitting_error, diff_constr,

                         nlp_solver_func);

        std::cout << "OK: Solver set up\n";


        /** Ident iteration loop =================================================================================== **/
        /* Print current sequence id (session/flight/seq */
        std::cout << "Running " << nIt << " iteration(s) on base parameters: " << kite_baseParams_filename << "\n";

        std::string kite_params_out_filepath;
        DM temp_x;
        DM temp_lam_g;
        DM temp_lam_x;

        int i_optimization_failed{0};

        for (int iIt = 0; iIt < nIt; ++iIt) {

            /** Select parameter IN file for current iteration **/
            if (iIt > 0) {
                /* There was an iteration before, use iterated param file. */
                kite_params_in_filepath = kite_params_out_filepath;
            }

            std::string kite_params_out_dir = seq_dir + "params_new/" + kite_baseParams_filename + "/";
            std::string afterItStr = "afterIt_" + std::to_string(iIt + 1);
            kite_params_out_filepath = kite_params_out_dir + afterItStr + ".yaml";

            if (kite_params_warmStart) {
                if (kite_utils::file_exists(kite_params_out_filepath)) {
                    continue;
                }
            }

            std::cout << "\nIteration " << iIt + 1 << "\n";
            std::cout << "Kite param file IN: " << kite_params_in_filepath << "\n";
            KiteProperties kite_props_in = kite_utils::LoadProperties(kite_params_in_filepath);

            /* Parameter bounds */
            /* (For each DATA_POINT) */
            DM REF_P;   // Parameter initial values
            DM LBP;     // Parameter lower bounds
            DM UBP;     // Parameter upper bounds

            std::list<Parameter> parameterList; // List to map parameter id and name to sort them by sensitivity later

            setup_optimizationParameters(kite_props_in, identManeuver,

                                         REF_P, LBP, UBP, parameterList);
            std::cout << "Parameter preparation OK\n";


            /** set default args */
            DMDict ARG;
            SX lbx, ubx;

            get_optimizationVarBounds(LBX, UBX, num_segments, poly_order,
                                      LBU, UBU,
                                      LBP, UBP,

                                      lbx, ubx);

            SX lbg = SX::zeros(diff_constr.size());
            SX ubg = SX::zeros(diff_constr.size());


            /* Feasible optimization variable guess  */
            DM feasible_state;

            DMDict feasible_guess;

            if (!temp_x.is_empty()) {
                /* Use found trajectory from last iteration as initial feasible guess for this iteration */
                std::cout << "Using initial guess from last iteration.\n";

                ARG["x0"] = temp_x;
                ARG["lam_x0"] = temp_lam_x;
                ARG["lam_g0"] = temp_lam_g;

            } else {
                /* Provide initial guess from integrator by flying at constant control input */
                std::cout << "Computing initial guess from integrator.\n";

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

                std::cout << "init_state: " << id_state0 << " init_control: " << init_control << "\n";

                feasible_guess = ps_solver.solve_trajectory(id_state0, init_control, true);

                DM feasible_traj = feasible_guess.at("x");
                ARG["x0"] = casadi::DM::vertcat({feasible_traj, REF_P});
                ARG["lam_x0"] = DM::vertcat({feasible_guess.at("lam_x"), DM::zeros(REF_P.size1())});
                ARG["lam_g0"] = feasible_guess.at("lam_g");
                std::cout << "Initial guess found.\n";
            }

            /* Set initial state
             * As the problem is formulated time-reverse, the position of the initial state
             * in the optimization variable vector is at the end of the states. */
            int idx_in = num_segments * poly_order * dimx;
            int idx_out = idx_in + dimx;
            lbx(Slice(idx_in, idx_out), 0) = id_state0;
            ubx(Slice(idx_in, idx_out), 0) = id_state0;

            /** Solve the optimization problem ===================================================================== **/
            std::cout << "\nSolving the optimization problem ...\n";

            ARG["lbx"] = lbx;
            ARG["ubx"] = ubx;
            ARG["lbg"] = lbg;
            ARG["ubg"] = ubg;

            DMDict res = nlp_solver_func(ARG);

            /* Check if optimization was infeasible. */
            Dict stats = nlp_solver_func.stats();
            std::string solve_status = static_cast<std::string>(stats["return_status"]);
            if (solve_status == "Invalid_Number_Detected") {
                ++i_optimization_failed;
                if (i_optimization_failed >= n_optimizationFailed_max) {
                    std::cout << "Invalid number detected!\n";
                    iIt = nIt;
                    continue;
                }
            }
            if (solve_status == "Infeasible_Problem_Detected") {
                ++i_optimization_failed;
                if (i_optimization_failed >= n_optimizationFailed_max) {
                    std::cout << "Infeasible problem detected. Stopping iterations on this sequence. \n";
                    iIt = nIt;
                    continue;
                }
            }
            if (solve_status == "Restoration_Failed") {
                ++i_optimization_failed;
                if (i_optimization_failed >= n_optimizationFailed_max) {
                    std::cout << "Restoration failed! Stopping iterations on this sequence. \n";
                    iIt = nIt;
                    continue;
                }
            }

            /* Save results for next iteration */
            temp_x = res.at("x");
            temp_lam_x = res.at("lam_x");
            temp_lam_g = res.at("lam_g");

            DM result = res.at("x");
            DM lam_x = res.at("lam_x");

            /** Output ============================================================================================= **/
            /** Write new parameters to file **/
            /* Extract optimized parameters from optimization result */
            DM new_params = result(Slice(result.size1() - varp.size1(), result.size1()));
            std::vector<double> new_params_vec = new_params.nonzeros();

            write_parameterFile(kite_params_in_filepath, kite_params_out_filepath,
                                identManeuver, new_params_vec);

            /** Write estimated trajectory to file **/
            std::string est_traj_filepath = kite_params_out_dir + afterItStr + "_estimatedTrajectory.txt";

            /* Time-reverse trajectory from optimization result, reshape, reverse time */
            DM est_traj_rev_vect = result(Slice(0, varx.size1()));
            DM est_traj_rev_woTime = DM::reshape(est_traj_rev_vect, dimx, DATA_POINTS);
            DM est_traj_woTime = flipDM(est_traj_rev_woTime, 2);

            /* Now: Rows: States, Columns: Time steps
             * Add original times as first row */
            DM est_traj_wTime = DM::vertcat({id_time,
                                             est_traj_woTime});

            write_trajectoryFile(est_traj_wTime, est_traj_filepath);

            /** Visualize parameters found within their bounds and write to file **/
            /* Get parameter sensitivities*/
            DM param_sens = lam_x(Slice(lam_x.size1() - varp.size1(), lam_x.size1()));

            /* Calculate fitting error (cost function) */
            DM fitting_error_evaluated = fitting_error_func(DMVector{est_traj_rev_woTime})[0];

            std::string paramInfo_filepath = kite_params_out_dir + afterItStr + "_paramInfo.txt";

            printWrite_parametersFound(parameterList, new_params_vec, LBP, UBP, param_sens, fitting_error_evaluated,
                                       paramInfo_filepath);

        } /** End of ident iteration loop -------------------------------------------------------------------------- **/

    } /** End of sequence loop ------------------------------------------------------------------------------------- **/

}
