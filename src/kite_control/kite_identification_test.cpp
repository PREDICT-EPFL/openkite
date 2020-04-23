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

    int id{0};
    const std::string groupName;
    const std::string name;
    double refValue{};
    double lowerBound{};
    double upperBound{};

    double foundValue{};
    double sensitivity{};

    /* With groupName */
    Parameter(
            std::string groupName_,
            std::string name_,
            const double &refValue_,
            const double &lowerBound_,
            const double &upperBound_,
            const bool absBounds = false) :

            groupName(std::move(groupName_)),
            name(std::move(name_)),
            refValue(refValue_),
            lowerBound((absBounds) ? lowerBound_ : (refValue_ + lowerBound_ * std::abs(refValue_))),
            upperBound((absBounds) ? upperBound_ : (refValue_ + upperBound_ * std::abs(refValue_))) {}

    Parameter(std::string name_) :

            name(std::move(name_)) {}

    Parameter() = default;


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


void readin_optResultsFile(const std::string &filepath,

                           DM &optResult_x, DM &optResult_lam_x, DM &optResult_lam_g) {

    std::ifstream optRes_file(filepath, std::ios::in);

    if (!optRes_file.fail()) {


        std::stringstream linestream;
        std::string line;

        std::string str;
        double number{0};

        /** First line **/
        getline(optRes_file, line);
        linestream << line;

        linestream >> str;
        int i{0};
        while (linestream >> number) {
            optResult_x(i) = number;
            ++i;
        }
        linestream.clear();

        /** Second line **/
        getline(optRes_file, line);
        linestream << line;

        linestream >> str;
        int j{0};
        while (linestream >> number) {
            optResult_lam_x(j) = number;
            ++j;
        }
        linestream.clear();

        /** Third line **/
        getline(optRes_file, line);
        linestream << line;

        linestream >> str;
        int k{0};
        while (linestream >> number) {
            optResult_lam_g(k) = number;
            ++k;
        }

    } else {
        std::cout << "Could not open optResult file.\n";
        optRes_file.clear();
    }
}


/* Set up */
void get_costMatrix(const kite_utils::IdentMode identMode, bool useEulerAngles,

                    DM &Q) {

    SX scalingMat;

    if (useEulerAngles) {

        /* Default */
        Q = SX::diag(SX({100, 100, 100,   // vx vy vz
                         200, 200, 200,   // wx wy wz
                         10, 10, 10,       // x y z
                         100, 100, 100})); //roll qy yaw

        if (identMode == kite_utils::IdentMode::LONGITUDINAL) {

            Q = SX::diag(SX({10, 10, 100,       // vx vy vz
                             10, 200, 10,      // wx wy wz
                             10, 10, 100,       // x y z
                             10, 200, 10}));     //roll qy yaw

        } else if (identMode == kite_utils::IdentMode::LATERAL) {

            Q = SX::diag(SX({10, 100, 10,   // vx vy vz
                             100, 10, 100,   // wx wy wz
                             10, 10, 10,   // x y z
                             100, 10, 100}));   //roll qy yaw

        } else if (identMode == kite_utils::IdentMode::YAW) {

            Q = SX::diag(SX({100, 50, 10,   // vx vy vz
                             100, 10, 100,   // wx wy wz
                             10, 10, 10,   // x y z
                             100, 10, 100}));   //roll qy yaw

        } else {

            Q = SX::diag(SX::ones(12) * 100);
        }

        /* Scaling */
        scalingMat = SX::diag(SX({1.0 / 20, 1.0 / 3, 1.0 / 5,   // vx vy vz
                                  1.0 / 3.14, 1.0 / 3.14, 1.0 / 3.14,   // wx wy wz
                                  1.0 / 30, 1.0 / 30, 1.0 / 30,   // x y z
                                  1.0 / 1.047, 1.0 / 1, 1.0 / 0.785}));   //roll qy yaw

    } else {
        /* Default */
        Q = SX::diag(SX({100, 100, 100,   // vx vy vz
                         200, 200, 200,   // wx wy wz
                         10, 10, 10,       // x y z
                         100, 100, 100, 100})); //qw qx qy qz

        if (identMode == kite_utils::IdentMode::LONGITUDINAL) {

            Q = SX::diag(SX({100, 10, 100,   // vx vy vz
                             10, 100, 10,   // wx wy wz
                             10, 10, 100,   // x y z
                             10, 10, 100, 10}));   //qw qx qy qz

        } else if (identMode == kite_utils::IdentMode::LATERAL) {

            Q = SX::diag(SX({10, 100, 10,   // vx vy vz
                             100, 10, 100,   // wx wy wz
                             10, 10, 10,   // x y z
                             10, 10, 10, 10}));   //qw qx qy qz

        } else if (identMode == kite_utils::IdentMode::YAW) {

            Q = SX::diag(SX({100, 50, 10,   // vx vy vz
                             100, 10, 100,   // wx wy wz
                             10, 10, 10,   // x y z
                             10, 10, 10, 10}));   //qw qx qy qz

        } else {

            Q = SX::diag(SX::ones(13) * 100);
        }

        /* Scaling */
        scalingMat = SX::diag(SX({1.0 / 20, 1.0 / 3, 1.0 / 5,   // vx vy vz
                                  1.0 / 3.14, 1.0 / 3.14, 1.0 / 3.14,   // wx wy wz
                                  1.0 / 30, 1.0 / 30, 1.0 / 30,   // x y z
                                  1.0 / 1, 1.0 / 1, 1.0 / 1, 1.0 / 1}));   //qw qx qy qz
    }

    Q = SX::mtimes(scalingMat, Q);

}

void get_kiteDynamics(KiteProperties &kite_props, const double &windFrom_deg, const double &windSpeed,
                      const kite_utils::IdentMode &identMode, const bool controlsIncludeWind,

                      KiteDynamics &kite, KiteDynamics &kite_int) {

    kite_props.Wind.WindFrom = windFrom_deg * M_PI / 180.0;
    kite_props.Wind.WindSpeed = windSpeed;

    AlgorithmProperties algo_props;
    algo_props.Integrator = CVODES;
    algo_props.sampling_time = 0.02;

    kite = KiteDynamics(kite_props, algo_props, identMode, controlsIncludeWind);
    kite_int = KiteDynamics(kite_props, algo_props, controlsIncludeWind); //integration model
}

void get_minimalKiteDynamics(KiteProperties &kite_props, const double &windFrom_deg, const double &windSpeed,
                             const kite_utils::IdentMode &identMode, const bool controlsIncludeWind,

                             MinimalKiteDynamics &kite, MinimalKiteDynamics &kite_int) {

    kite_props.Wind.WindFrom = windFrom_deg * M_PI / 180.0;
    kite_props.Wind.WindSpeed = windSpeed;

    AlgorithmProperties algo_props;
    algo_props.Integrator = CVODES;
    algo_props.sampling_time = 0.02;

    kite = MinimalKiteDynamics(kite_props, algo_props, identMode, controlsIncludeWind);
    kite_int = MinimalKiteDynamics(kite_props, algo_props, controlsIncludeWind); //integration model
}

void substitute_quat2euler(SX &state) {

    SX qw = state(9);
    SX qx = state(10);
    SX qy = state(11);
    SX qz = state(12);

    /* MATLAB quat2eul function, Copyright 2014-2017 The MathWorks, Inc. */

    // Not differentiable -> Keep qy instead
//    SX aSinInput = -2 * (qx * qz - qw * qy);
//    aSinInput = SX::fmin(aSinInput, 1);
//    aSinInput = SX::fmax(aSinInput, -1);
//    SX pitch = asin(aSinInput);

    SX yaw = atan2(2 * (qx * qy + qw * qz), qw * qw + qx * qx - qy * qy - qz * qz);
    SX roll = atan2(2 * (qy * qz + qw * qx), qw * qw - qx * qx - qy * qy + qz * qz);

    state = SX::vertcat({state(Slice(0, 9)), roll, qy, yaw});
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

void setup_optimizationParameters(const KiteProperties &kiteProps,
                                  const kite_utils::IdentMode &identMode, const bool setBounds,

                                  DM &REF_P, DM &LBP, DM &UBP, std::list<Parameter> &paramList) {

    /* Define optimization parameters and their bounds.
     * If setBounds = false, the relative or absolute bounds set to each parameter will not be written at the end.
     * They only refer to the first time this function is called and therefore the base kite parameters */

    if (identMode == kite_utils::IdentMode::LONGITUDINAL) {

        /** All at once, allowing negative Cm0 (due to high wing setting angle) **/
        const double relBound = 0.7;

        paramList.emplace_back("aero", "CD0", kiteProps.Aerodynamics.CD0, -relBound, relBound);

        /* AOA */
        paramList.emplace_back("aero_aoa", "CL0", kiteProps.Aerodynamics.CL0, -relBound, relBound);
        paramList.emplace_back("aero_aoa", "CLa", kiteProps.Aerodynamics.CLa, -relBound, relBound);

        paramList.emplace_back("aero_aoa", "Cm0", kiteProps.Aerodynamics.Cm0, -0.1, 0.03, true);
        paramList.emplace_back("aero_aoa", "Cma", kiteProps.Aerodynamics.Cma, -relBound, relBound);


        /* Pitchrate */
        paramList.emplace_back("aero_rate_pitch", "CLq", kiteProps.Aerodynamics.CLq, -relBound, relBound);
        paramList.emplace_back("aero_rate_pitch", "Cmq", kiteProps.Aerodynamics.Cmq, -relBound, relBound);


        /* Elevator */
        paramList.emplace_back("aero_ctrl_elev", "CLde", kiteProps.Aerodynamics.CLde, -0.3, 0.3);
        paramList.emplace_back("aero_ctrl_elev", "Cmde", kiteProps.Aerodynamics.Cmde, -0.3, 0.3);
        // 9 longitudinal parameters

    } else if (identMode == kite_utils::IdentMode::LATERAL) {

        double relUpBound = 0.7;
        double relLwBound = 0.9;

        /* Sideslip */
        paramList.emplace_back("aero_ss", "CYb", kiteProps.Aerodynamics.CYb, -relLwBound, relUpBound);

        paramList.emplace_back("aero_ss", "Clb", kiteProps.Aerodynamics.Clb, -relLwBound, relUpBound);

        paramList.emplace_back("aero_ss", "Cn0", kiteProps.Aerodynamics.Cn0, -relLwBound, relUpBound);
        paramList.emplace_back("aero_ss", "Cnb", kiteProps.Aerodynamics.Cnb, -relLwBound, relUpBound);


        /* Rollrate */
        paramList.emplace_back("aero_rate_roll", "CYp", kiteProps.Aerodynamics.CYp, -relLwBound, relUpBound);
        paramList.emplace_back("aero_rate_roll", "Clp", kiteProps.Aerodynamics.Clp, -relLwBound, relUpBound);
        paramList.emplace_back("aero_rate_roll", "Cnp", kiteProps.Aerodynamics.Cnp, -relLwBound, relUpBound);

        /* Yawrate */
        paramList.emplace_back("aero_rate_yaw", "CYr", kiteProps.Aerodynamics.CYr, -relLwBound, relUpBound);
        paramList.emplace_back("aero_rate_yaw", "Clr", kiteProps.Aerodynamics.Clr, -relLwBound, relUpBound);
        paramList.emplace_back("aero_rate_yaw", "Cnr", kiteProps.Aerodynamics.Cnr, -relLwBound, relUpBound);


        /* Ailerons */
        paramList.emplace_back("aero_ctrl_ail", "Clda", kiteProps.Aerodynamics.Clda, -relLwBound, relUpBound);
        paramList.emplace_back("aero_ctrl_ail", "Cnda", kiteProps.Aerodynamics.Cnda, -relLwBound, relUpBound);
        // 12 lateral/aileron parameters


    } else if (identMode == kite_utils::IdentMode::YAW) {

        double relUpBound = 0.7;
        double relLwBound = 0.9;

        /* Sideslip */
        paramList.emplace_back("aero_ss", "CYb", kiteProps.Aerodynamics.CYb, -relLwBound, relUpBound);

        paramList.emplace_back("aero_ss", "Clb", kiteProps.Aerodynamics.Clb, -relLwBound, relUpBound);

        paramList.emplace_back("aero_ss", "Cn0", kiteProps.Aerodynamics.Cn0, 0, 0); //fixed
        paramList.emplace_back("aero_ss", "Cnb", kiteProps.Aerodynamics.Cnb, -relLwBound, relUpBound);


        /* Rollrate */
        paramList.emplace_back("aero_rate_roll", "CYp", kiteProps.Aerodynamics.CYp, -relLwBound, relUpBound);
        paramList.emplace_back("aero_rate_roll", "Clp", kiteProps.Aerodynamics.Clp, -0, 0);
        paramList.emplace_back("aero_rate_roll", "Cnp", kiteProps.Aerodynamics.Cnp, -0, 0);

        /* Yawrate */
        paramList.emplace_back("aero_rate_yaw", "CYr", kiteProps.Aerodynamics.CYr, -relLwBound, relUpBound);
        paramList.emplace_back("aero_rate_yaw", "Clr", kiteProps.Aerodynamics.Clr, -0, 0);
        paramList.emplace_back("aero_rate_yaw", "Cnr", kiteProps.Aerodynamics.Cnr, -0, 0);


        /* Rudder */
        paramList.emplace_back("aero_ctrl_rud", "CYdr", kiteProps.Aerodynamics.CYdr, -relLwBound, relUpBound);
        paramList.emplace_back("aero_ctrl_rud", "Cldr", kiteProps.Aerodynamics.Cldr, -relLwBound, relUpBound);
        paramList.emplace_back("aero_ctrl_rud", "Cndr", kiteProps.Aerodynamics.Cndr, -relLwBound, relUpBound);
        // 13 lateral/rudder parameters

    } else if (identMode == kite_utils::IdentMode::COMPLETE) {

        double relBound = 0.7;

        paramList.emplace_back("aero", "CD0", kiteProps.Aerodynamics.CD0, -relBound, relBound);

        /* AOA */
        paramList.emplace_back("aero_aoa", "CL0", kiteProps.Aerodynamics.CL0, -relBound, relBound);
        paramList.emplace_back("aero_aoa", "CLa", kiteProps.Aerodynamics.CLa, -relBound, relBound);

        paramList.emplace_back("aero_aoa", "Cm0", kiteProps.Aerodynamics.Cm0, -0.03, 0.03, true);
        paramList.emplace_back("aero_aoa", "Cma", kiteProps.Aerodynamics.Cma, -relBound, relBound);

        /* Sideslip */
        paramList.emplace_back("aero_ss", "CYb", kiteProps.Aerodynamics.CYb, -relBound, relBound);

        paramList.emplace_back("aero_ss", "Clb", kiteProps.Aerodynamics.Clb, -relBound, relBound);

        paramList.emplace_back("aero_ss", "Cn0", kiteProps.Aerodynamics.Cn0, -0.2, 0.2, true);
        paramList.emplace_back("aero_ss", "Cnb", kiteProps.Aerodynamics.Cnb, -relBound, relBound);


        /* Pitchrate */
        paramList.emplace_back("aero_rate_pitch", "CLq", kiteProps.Aerodynamics.CLq, -relBound, relBound);
        paramList.emplace_back("aero_rate_pitch", "Cmq", kiteProps.Aerodynamics.Cmq, -relBound, relBound);

        /* Rollrate */
        paramList.emplace_back("aero_rate_roll", "CYp", kiteProps.Aerodynamics.CYp, -relBound, relBound);
        paramList.emplace_back("aero_rate_roll", "Clp", kiteProps.Aerodynamics.Clp, -relBound, relBound);
        paramList.emplace_back("aero_rate_roll", "Cnp", kiteProps.Aerodynamics.Cnp, -relBound, relBound);

        /* Yawrate */
        paramList.emplace_back("aero_rate_yaw", "CYr", kiteProps.Aerodynamics.CYr, -relBound, relBound);
        paramList.emplace_back("aero_rate_yaw", "Clr", kiteProps.Aerodynamics.Clr, -relBound, relBound);
        paramList.emplace_back("aero_rate_yaw", "Cnr", kiteProps.Aerodynamics.Cnr, -relBound, relBound);


        /* Elevator */
        paramList.emplace_back("aero_ctrl_elev", "CLde", kiteProps.Aerodynamics.CLde, -relBound, relBound);
        paramList.emplace_back("aero_ctrl_elev", "Cmde", kiteProps.Aerodynamics.Cmde, -relBound, relBound);

        /* Ailerons */
        paramList.emplace_back("aero_ctrl_ail", "Clda", kiteProps.Aerodynamics.Clda, -relBound, relBound);
        paramList.emplace_back("aero_ctrl_ail", "Cnda", kiteProps.Aerodynamics.Cnda, -relBound, relBound);

        /* Rudder */
        paramList.emplace_back("aero_ctrl_rud", "CYdr", kiteProps.Aerodynamics.CYdr, -relBound, relBound);
        paramList.emplace_back("aero_ctrl_rud", "Cldr", kiteProps.Aerodynamics.Cldr, -relBound, relBound);
        paramList.emplace_back("aero_ctrl_rud", "Cndr", kiteProps.Aerodynamics.Cndr, -relBound, relBound);

        // 25 complete parameters
    }


    if (setBounds) {
        /* Setting bounds for the first time */
        LBP = REF_P = UBP = DM::zeros(paramList.size());

        /* Numerize parameters and fill bounds */
        int i = 0;
        for (Parameter &param : paramList) {
            param.id = i;
            UBP(i) = param.upperBound;
            REF_P(i) = param.refValue;
            LBP(i) = param.lowerBound;

            ++i;
        }
        std::cout << "Base parameters set.\n";

    } else {
        /* Reuse bounds from previous iteration (do not overwrite them).  */
        /* Numerize parameters and fill bounds */
        int i = 0;
        for (Parameter &param : paramList) {
            param.id = i;
            param.upperBound = UBP(i).nonzeros()[0];
            REF_P(i) = param.refValue;
            param.lowerBound = LBP(i).nonzeros()[0];

            ++i;
        }
    }

//    std::cout << "Upper Bound: " << UBP << "\n"
//              << "Ref value  : " << REF_P << "\n"
//              << "Lower bound: " << LBP << "\n";
}

void setup_minimalOptimizationParameters(const KiteProperties &kiteProps,
                                         const kite_utils::IdentMode &identMode, const bool setBounds,

                                         DM &REF_P, DM &LBP, DM &UBP, std::list<Parameter> &paramList) {

    /* Define optimization parameters and their bounds.
     * If setBounds = false, the relative or absolute bounds set to each parameter will not be written at the end.
     * They only refer to the first time this function is called and therefore the base kite parameters */

    if (identMode == kite_utils::IdentMode::LONGITUDINAL) {

        /** All at once, allowing negative Cm0 (due to high wing setting angle) **/
        const double relBound = 0.7;

        paramList.emplace_back("aero", "CD0", kiteProps.Aerodynamics.CD0, -relBound, relBound);

        /* AOA */
        paramList.emplace_back("aero_aoa", "CL0", kiteProps.Aerodynamics.CL0, -relBound, relBound);
        paramList.emplace_back("aero_aoa", "CLa", kiteProps.Aerodynamics.CLa, -relBound, relBound);

        paramList.emplace_back("aero_aoa", "Cm0", kiteProps.Aerodynamics.Cm0, -0.03, 0.03, true);
        paramList.emplace_back("aero_aoa", "Cma", kiteProps.Aerodynamics.Cma, -relBound, relBound);


        /* Pitchrate */
//        paramList.emplace_back("aero_rate_pitch", "CLq", kiteProps.Aerodynamics.CLq, -relBound, relBound);
//        paramList.emplace_back("aero_rate_pitch", "Cmq", kiteProps.Aerodynamics.Cmq, -relBound, relBound);


        /* Elevator */
//        paramList.emplace_back("aero_ctrl_elev", "CLde", kiteProps.Aerodynamics.CLde, -relBound, relBound);
        paramList.emplace_back("aero_ctrl_elev", "Cmde", kiteProps.Aerodynamics.Cmde, -relBound, relBound);
        // 6 longitudinal parameters

    } else if (identMode == kite_utils::IdentMode::LATERAL) {

        double relBound = 0.7;

        /* Sideslip */
        paramList.emplace_back("aero_ss", "CYb", kiteProps.Aerodynamics.CYb, -relBound, relBound);

        paramList.emplace_back("aero_ss", "Clb", kiteProps.Aerodynamics.Clb, -relBound, relBound);

        paramList.emplace_back("aero_ss", "Cnb", kiteProps.Aerodynamics.Cnb, -relBound, relBound);


        /* Rollrate */
        //paramList.emplace_back("aero_rate_roll", "CYp", kiteProps.Aerodynamics.CYp, -relBound, relBound);
        paramList.emplace_back("aero_rate_roll", "Clp", kiteProps.Aerodynamics.Clp, -relBound, relBound);
//        paramList.emplace_back("aero_rate_roll", "Cnp", kiteProps.Aerodynamics.Cnp, -relBound, relBound);

        /* Yawrate */
//        paramList.emplace_back("aero_rate_yaw", "CYr", kiteProps.Aerodynamics.CYr, -relBound, relBound);
//        paramList.emplace_back("aero_rate_yaw", "Clr", kiteProps.Aerodynamics.Clr, -relBound, relBound);
        paramList.emplace_back("aero_rate_yaw", "Cnr", kiteProps.Aerodynamics.Cnr, -relBound, relBound);


        /* Ailerons */
        paramList.emplace_back("aero_ctrl_ail", "Clda", kiteProps.Aerodynamics.Clda, -relBound, relBound);
        //paramList.emplace_back("aero_ctrl_ail", "Cnda", kiteProps.Aerodynamics.Cnda, -relBound, relBound);
        // 6 lateral/aileron parameters


    } else if (identMode == kite_utils::IdentMode::YAW) {

        double relBound = 0.7;

        /* Sideslip */
        paramList.emplace_back("aero_ss", "CYb", kiteProps.Aerodynamics.CYb, -relBound, relBound);

        paramList.emplace_back("aero_ss", "Clb", kiteProps.Aerodynamics.Clb, -relBound, relBound);

        paramList.emplace_back("aero_ss", "Cnb", kiteProps.Aerodynamics.Cnb, -relBound, relBound);


        /* Rollrate */
        //paramList.emplace_back("aero_rate_roll", "CYp", kiteProps.Aerodynamics.CYp, -relBound, relBound);
        paramList.emplace_back("aero_rate_roll", "Clp", kiteProps.Aerodynamics.Clp, -relBound, relBound);
//        paramList.emplace_back("aero_rate_roll", "Cnp", kiteProps.Aerodynamics.Cnp, -relBound, relBound);

        /* Yawrate */
//        paramList.emplace_back("aero_rate_yaw", "CYr", kiteProps.Aerodynamics.CYr, -relBound, relBound);
//        paramList.emplace_back("aero_rate_yaw", "Clr", kiteProps.Aerodynamics.Clr, -relBound, relBound);
        paramList.emplace_back("aero_rate_yaw", "Cnr", kiteProps.Aerodynamics.Cnr, -relBound, relBound);


        /* Rudder */
//        paramList.emplace_back("aero_ctrl_rud", "CYdr", kiteProps.Aerodynamics.CYdr, -relBound, relBound);
//        paramList.emplace_back("aero_ctrl_rud", "Cldr", kiteProps.Aerodynamics.Cldr, -relBound, relBound);
        paramList.emplace_back("aero_ctrl_rud", "Cndr", kiteProps.Aerodynamics.Cndr, -relBound, relBound);
        // 6 lateral/rudder parameters

    }

//    else if (identMode == kite_utils::IdentMode::COMPLETE) {
//
//        paramList.emplace_back("aero", "CD0", kiteProps.Aerodynamics.CD0, -0.5, 0.7);
//
//        /* AOA */
//        paramList.emplace_back("aero_aoa", "CL0", kiteProps.Aerodynamics.CL0, -0.1, 0.5);
//        paramList.emplace_back("aero_aoa", "CLa", kiteProps.Aerodynamics.CLa, -0.05, 0.1);
//
//        paramList.emplace_back("aero_aoa", "Cm0", kiteProps.Aerodynamics.Cm0, -0.5, 0.5);
//        paramList.emplace_back("aero_aoa", "Cma", kiteProps.Aerodynamics.Cma, -0.1, 0.3);
//
//        /* Sideslip */
//        paramList.emplace_back("aero_ss", "CYb", kiteProps.Aerodynamics.CYb, -0.5, 0.5);
//
//        paramList.emplace_back("aero_ss", "Clb", kiteProps.Aerodynamics.Clb, -0.5, 0.5);
//
//        paramList.emplace_back("aero_ss", "Cnb", kiteProps.Aerodynamics.Cnb, -0.5, 0.5);
//
//
//        /* Pitchrate */
//        paramList.emplace_back("aero_rate_pitch", "CLq", kiteProps.Aerodynamics.CLq, -0.5, 0.5);
//        paramList.emplace_back("aero_rate_pitch", "Cmq", kiteProps.Aerodynamics.Cmq, -0.3, 0.3);
//
//        /* Rollrate */
//        paramList.emplace_back("aero_rate_roll", "CYp", kiteProps.Aerodynamics.CYp, -0.5, 0.5);
//        paramList.emplace_back("aero_rate_roll", "Clp", kiteProps.Aerodynamics.Clp, -0.5, 0.5);
//        paramList.emplace_back("aero_rate_roll", "Cnp", kiteProps.Aerodynamics.Cnp, -0.3, 1.0);
//
//        /* Yawrate */
//        paramList.emplace_back("aero_rate_yaw", "CYr", kiteProps.Aerodynamics.CYr, -0.3, 0.3);
//        paramList.emplace_back("aero_rate_yaw", "Clr", kiteProps.Aerodynamics.Clr, -0.5, 0.5);
//        paramList.emplace_back("aero_rate_yaw", "Cnr", kiteProps.Aerodynamics.Cnr, -0.3, 1.0);
//
//
//        /* Elevator */
//        paramList.emplace_back("aero_ctrl_elev", "CLde", kiteProps.Aerodynamics.CLde, -0.5, 0.5);
//        paramList.emplace_back("aero_ctrl_elev", "Cmde", kiteProps.Aerodynamics.Cmde, -0.5, 0.5);
//
//        /* Ailerons */
//        paramList.emplace_back("aero_ctrl_ail", "Clda", kiteProps.Aerodynamics.Clda, -0.5, 0.5);
//        paramList.emplace_back("aero_ctrl_ail", "Cnda", kiteProps.Aerodynamics.Cnda, -0.5, 1.5);
//        // 21 complete parameters
//
////        /* Rudder */
////        paramList.emplace_back("aero_ctrl_rud", "CYdr", kiteProps.Aerodynamics.CYdr, -0.5, 0.5);
////        paramList.emplace_back("aero_ctrl_rud", "Cldr", kiteProps.Aerodynamics.Cldr, -0.5, 0.5);
////        paramList.emplace_back("aero_ctrl_rud", "Cndr", kiteProps.Aerodynamics.Cndr, -0.5, 0.5);
//    }


    if (setBounds) {
        /* Setting bounds for the first time */
        LBP = REF_P = UBP = DM::zeros(paramList.size());

        /* Numerize parameters and fill bounds */
        int i = 0;
        for (Parameter &param : paramList) {
            param.id = i;
            UBP(i) = param.upperBound;
            REF_P(i) = param.refValue;
            LBP(i) = param.lowerBound;

            ++i;
        }
        std::cout << "Base parameters set.\n";

    } else {
        /* Reuse bounds from previous iteration (do not overwrite them).  */
        /* Numerize parameters and fill bounds */
        int i = 0;
        for (Parameter &param : paramList) {
            param.id = i;
            param.upperBound = UBP(i).nonzeros()[0];
            REF_P(i) = param.refValue;
            param.lowerBound = LBP(i).nonzeros()[0];

            ++i;
        }
    }

//    std::cout << "Upper Bound: " << UBP << "\n"
//              << "Ref value  : " << REF_P << "\n"
//              << "Lower bound: " << LBP << "\n";
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
void write_parameterFile(const std::string &kite_params_in_filepath, const std::list<Parameter> &paramList,
                         const std::string &kite_params_out_filepath) {

    /** update parameter file */
    YAML::Node kite_params_yaml = YAML::LoadFile(kite_params_in_filepath);

    for (auto &param : paramList) {
        kite_params_yaml[param.groupName][param.name] = param.foundValue;
    }

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

void write_optResultsFile(const DM &optResult_x, const DM &optResult_lam_x, const DM &optResult_lam_g,
                          const std::string &filepath) {

    std::ofstream optRes_file(filepath, std::ios::out);
    optRes_file.precision(15);

    optRes_file << "optResult_x: ";
    for (int i = 0; i < optResult_x.size1() - 1; i++) {
        optRes_file << optResult_x(i) << " ";
    }
    optRes_file << optResult_x(optResult_x.size1() - 1) << "\n";

    optRes_file << "optResult_lam_x: ";
    for (int i = 0; i < optResult_lam_x.size1() - 1; i++) {
        optRes_file << optResult_lam_x(i) << " ";
    }
    optRes_file << optResult_lam_x(optResult_lam_x.size1() - 1) << "\n";

    optRes_file << "optResult_lam_g: ";
    for (int i = 0; i < optResult_lam_g.size1() - 1; i++) {
        optRes_file << optResult_lam_g(i) << " ";
    }
    optRes_file << optResult_lam_g(optResult_lam_g.size1() - 1);

    optRes_file.close();
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

void printWrite_parametersFound(std::list<Parameter> &paramList, const double &fitting_error_evaluated,
                                const std::string &solve_status,
                                const std::string &filepath) {

    /* Sort parameter list by sensitivity */
    paramList.sort();

    double sensitivity_max = paramList.begin()->sensitivity;

    /* Header */
    std::cout << "\n";
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
    for (Parameter &param : paramList) {
        int i = param.id;

        print_singleParameterOutput(i, param.name,
                                    param.foundValue,
                                    param.lowerBound,
                                    param.upperBound,
                                    param.sensitivity, sensitivity_max);
    }
    std::cout << "\nFitting error = " << fitting_error_evaluated << "\n";


    /** Save paramInfo file **/

    std::ofstream param_diag_file(filepath, std::ios::out);

    param_diag_file << "ParamId" << " "
                    << "GroupName" << " "
                    << "ParamName" << " "
                    << "ParamValue" << " "
                    << "LwBound" << " "
                    << "UpBound" << " "
                    << "Sensitivity" << " "
                    << "SolveSucceeded" << " "
                    << "FittingError" << "\n";

    bool solveSucceeded = (solve_status == "Solve_Succeeded");

    for (Parameter &param : paramList) {
        int i = param.id;

        param_diag_file << std::fixed << std::setprecision(15);

        param_diag_file << i << " "
                        << param.groupName << " "
                        << param.name << " "
                        << param.foundValue << " "
                        << param.lowerBound << " "
                        << param.upperBound << " "
                        << param.sensitivity << " "
                        << solveSucceeded << " "
                        << fitting_error_evaluated << "\n";

    }
    param_diag_file.close();

}


int main() {

    /// 5 lines to adapt by the user =============================================================================== ///

    std::string flightDataDir = "/home/johannes/identification/processed_flight_data/";

    /// State and Control dimensions ///
    const int dimx = 13;  // v(3) w(3) r(3) q(4)

    const bool controlsIncludeWind = false;
    const int dimu = 4;   // T elev rud ail + windFrom windSpeed

    /// 1. Identification mode ///
    /// For minimal model, also comment... #minimalmodel ///

//    // roll / lateral
//    const kite_utils::IdentMode identMode = kite_utils::IdentMode::LATERAL;
//    const int dimp = 12; // full model
////    const int dimp = 6; // minimal model
//    const bool useMinimalModel = false;
//
//    const bool useEulerAnglesForCostFunction = true;
//    const int DATA_POINTS = 127;
//    const int poly_order = 3;
//    const int num_segments = 42;


//    // pitch / longitudinal
//    const kite_utils::IdentMode identMode = kite_utils::IdentMode::LONGITUDINAL;
//    const int dimp = 9; // full model
////    const int dimp = 7; // minimal model
//    const bool useMinimalModel = false;
//
//    const bool useEulerAnglesForCostFunction = false;
//    const int DATA_POINTS = 148;
//    const int poly_order = 3;
//    const int num_segments = 49;

//    const int DATA_POINTS = 1261;
//    const int num_segments = 420;

    // yaw
    const kite_utils::IdentMode identMode = kite_utils::IdentMode::YAW;
    const int dimp = 13; // full model
//    const int dimp = 6; // minimal model
    const bool useMinimalModel = false;

    const bool useEulerAnglesForCostFunction = true;
    const int DATA_POINTS = 106;
    const int poly_order = 3;
    const int num_segments = 35;


//    // complete
//    const kite_utils::IdentMode identMode = kite_utils::IdentMode::COMPLETE;
//    const int dimp = 25; // full model
////    const int dimp = 6; // minimal model
//    const bool useMinimalModel = false;
//
//    const bool useEulerAnglesForCostFunction = false;
//    const int DATA_POINTS = 1261;
//    const int poly_order = 3;
//    const int num_segments = 420;


    const int n_optimizationFailed_max = 3;

    /// ============================================================================================================ ///
    std::string identManeuver;
    if (identMode == kite_utils::IdentMode::LONGITUDINAL) {
        identManeuver = "identPitch";
    } else if (identMode == kite_utils::IdentMode::LATERAL) {
        identManeuver = "identRoll";
    } else if (identMode == kite_utils::IdentMode::YAW) {
        identManeuver = "identYaw";
    } else if (identMode == kite_utils::IdentMode::COMPLETE) {
        identManeuver = "manual_id";
    }

    std::cout << "\nRunning " << identManeuver << " identification of " << dimp << " parameters with "
              << DATA_POINTS << " data points.\n";

    /* Load identification schedule */
    std::vector<FlightMetaData> flights = readIn_identSchedule(flightDataDir, identManeuver);

    /* Cost matrix */
    DM Q;
    get_costMatrix(identMode, useEulerAnglesForCostFunction,
                   Q);

    /* State bounds (constant over all sequences and iterations) */
    DM UBX = DM::vertcat({45, 10, 10,
                          1 * M_PI, 1 * M_PI, 1 * M_PI,
                          300, 300, 0,
                          1.05, 1.05, 1.05, 1.05});

    DM LBX = DM::vertcat({2.0, -10, -10,
                          -1 * M_PI, -1 * M_PI, -1 * M_PI,
                          -300, -300, -300,
                          -1.05, -1.05, -1.05, -1.05});     // for infinity use -DM::inf(1) and DM::inf(1)

    /** Collocation **/
    Chebyshev<SX, poly_order, num_segments, dimx, dimu, dimp> spectral;

    /* Vectors, length: states * collocation points */
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

            if (useEulerAnglesForCostFunction) {
                /* Convert quaternion to euler angles for cost calculation */
                substitute_quat2euler(measured_state);
                substitute_quat2euler(state_to_optimize);
            }
            /* Error and sum up */
            SX error = measured_state - state_to_optimize;
            fitting_error += static_cast<double>(1.0 / DATA_POINTS) * SX::sum1(SX::mtimes(Q, pow(error, 2)));

        }
        Function fitting_error_func = Function("fitting_error", {varx_}, {fitting_error});


        /** Kite Dynamics based on baseParams and current sequence's wind conditions **/
        /* Get KiteDynamics based on constant (during optimization) parameters.
         * Parameters to be optimized parameters are read in for each iteration. */
        std::string kite_params_in_filepath = kite_baseParams_dir + kite_baseParams_filename + ".yaml";
        std::cout << "Base params: " << kite_params_in_filepath << "\n";

        KiteProperties baseKiteParams;
        if (useMinimalModel)
            baseKiteParams = kite_utils::LoadMinimalProperties(kite_params_in_filepath);
        else
            baseKiteParams = kite_utils::LoadProperties(kite_params_in_filepath);

        /// *************** comment out appropriate model #minimalmodel ************** ///
        KiteDynamics kite, kite_int;
        get_kiteDynamics(baseKiteParams, windFrom_deg, windSpeed, identMode, controlsIncludeWind,
                         kite, kite_int);
//        MinimalKiteDynamics kite, kite_int;
//        get_minimalKiteDynamics(baseKiteParams, windFrom_deg, windSpeed, identMode, controlsIncludeWind,
//                                kite, kite_int);

        Function DynamicsFunc = kite.getNumericDynamics();
        //SX X = kite.getSymbolicState();
        //SX U = kite.getSymbolicControl();
        SX P = kite.getSymbolicParameters();

        Function ode = kite_int.getNumericDynamics();

        SX diff_constr = spectral.CollocateDynamics(DynamicsFunc, 0, tf);
        diff_constr = diff_constr(Slice(0, num_segments * poly_order * dimx));

        /** Set base parameters **/
        /* Parameter bounds */
        /* (For each DATA_POINT) */
        DM REF_P;   // Parameter ref values
        DM LBP;     // Parameter lower bounds
        DM UBP;     // Parameter upper bounds
        std::list<Parameter> paramList; // List to map parameter id and name to sort them by sensitivity later

        if (useMinimalModel)
            setup_minimalOptimizationParameters(baseKiteParams, identMode, true,
                                                REF_P, LBP, UBP, paramList);
        else
            setup_optimizationParameters(baseKiteParams, identMode, true,
                                         REF_P, LBP, UBP, paramList);

        /** Create solver **/
        Function nlp_solver_func;
        create_nlpSolver(opt_var, fitting_error, diff_constr,

                         nlp_solver_func);

        std::cout << "OK: Solver set up\n";


        /** Ident iteration loop =================================================================================== **/
        /* Print current sequence id (session/flight/seq */
        std::cout << "Running " << nIt << " iteration(s) on base parameters: " << kite_baseParams_filename << "\n";

        std::string kite_params_out_filepath;
        bool status_optimalSolutionFound{false};
        DM optResult_x;
        DM optResult_lam_x;
        DM optResult_lam_g;

        int i_optimization_failed{0};

        for (int iIt = 0; iIt < nIt; ++iIt) {

            /** Select parameter IN file for current iteration **/
            if (iIt > 0) {
                /* There was an iteration before, use iterated param file. */
                kite_params_in_filepath = kite_params_out_filepath;
            }

            /* Parameter out filepath after this iteration */
            std::string kite_params_out_dir = seq_dir + "params_new/" + kite_baseParams_filename + "/";
            std::string afterItStr = "afterIt_" + std::to_string(iIt + 1);
            kite_params_out_filepath = kite_params_out_dir + afterItStr + ".yaml";

            std::string optResults0_filepath = kite_params_out_dir + "optResults_0.txt";

            if (kite_params_warmStart && kite_utils::file_exists(kite_params_out_filepath)) {
                /* Skip this loop, as it has already been done. */
                std::cout << "\nSkipping parameter iteration " << iIt + 1 << " (already done).\n";
                continue;
            } /* Start iterated optimization from here. */

            std::cout << "\nParameter iteration " << iIt + 1 << "\n";

            std::cout << "Kite param file IN: " << kite_params_in_filepath << "\n";

            KiteProperties kite_props_in;
            if (useMinimalModel)
                kite_props_in = kite_utils::LoadMinimalProperties(kite_params_in_filepath);
            else
                kite_props_in = kite_utils::LoadProperties(kite_params_in_filepath);

            /* Update parameter ref values (result from last iteration) */
            paramList.clear();

            if (useMinimalModel)
                setup_minimalOptimizationParameters(kite_props_in, identMode, true,
                                                    REF_P, LBP, UBP, paramList);
            else
                setup_optimizationParameters(kite_props_in, identMode, true,
                                             REF_P, LBP, UBP, paramList);
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

            if (optResult_x.is_empty()) {
                if (kite_utils::file_exists(optResults0_filepath)) {
                    /* Guess is available from file */
                    std::cout << "Using initial guess from previous optimization (file).\n";

                    optResult_x = optResult_lam_x = DM::zeros(opt_var.size1());
                    optResult_lam_g = DM::zeros(diff_constr.size1());

                    readin_optResultsFile(optResults0_filepath,
                                          optResult_x, optResult_lam_x, optResult_lam_g);

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

                    DM init_control;
                    if (controlsIncludeWind)
                        init_control = DM({0.1, 0.0, 0.0, 0.0, 180.0, 0.5});
                    else
                        init_control = DM({0.1, 0.0, 0.0, 0.0});
                    //init_control = casadi::DM::repmat(init_control, (num_segments * poly_order + 1), 1);
                    init_control = DM::vec(id_controls_rev_woTime);
                    std::cout << "init_state: " << id_state0 << ",\n init control: " << init_control << "\n";

                    feasible_guess = ps_solver.solve_trajectory(id_state0, init_control, true);

                    DM feasible_traj = feasible_guess.at("x");
                    optResult_x = casadi::DM::vertcat({feasible_traj, REF_P});
                    optResult_lam_x = DM::vertcat({feasible_guess.at("lam_x"), DM::zeros(REF_P.size1())});
                    optResult_lam_g = feasible_guess.at("lam_g");
                    std::cout << "Initial guess found.\n";

                    /** Write estimated trajectory to file **/
                    /* Time-reverse trajectory from optimization result, reshape, reverse time */
                    DM est_traj_rev_vect = optResult_x(Slice(0, varx.size1()));
                    DM est_traj_rev_woTime = DM::reshape(est_traj_rev_vect, dimx, DATA_POINTS);
                    DM est_traj_woTime = flipDM(est_traj_rev_woTime, 2);

                    /* Now: Rows: States, Columns: Time steps
                     * Add original times as first row */
                    DM est_traj_wTime = DM::vertcat({id_time, est_traj_woTime});
                    std::string afterItStr = "afterIt_" + std::to_string(0);
                    std::string est_traj_filepath = kite_params_out_dir + afterItStr + "_estimatedTrajectory.txt";
                    write_trajectoryFile(est_traj_wTime, est_traj_filepath);
                }
            }

            /* Use trajectory from previous iteration as initial feasible guess for this iteration */
            ARG["x0"] = optResult_x;
            ARG["lam_x0"] = optResult_lam_x;
            ARG["lam_g0"] = optResult_lam_g;

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

            /* Check if optimization failed. */
            Dict stats = nlp_solver_func.stats();
            std::string solve_status = static_cast<std::string>(stats["return_status"]);

            if (solve_status == "Invalid_Number_Detected" ||
                solve_status == "Infeasible_Problem_Detected" ||
                solve_status == "Restoration_Failed") { ++i_optimization_failed; }

            /* Save results for next iteration */
            optResult_x = res.at("x");
            optResult_lam_x = res.at("lam_x");
            optResult_lam_g = res.at("lam_g");

            /** Output ============================================================================================= **/
            /* Get found parameter values and sensitivities from optimization result and write to paramList */
            DM new_params = optResult_x(Slice(optResult_x.size1() - varp.size1(), optResult_x.size1()));
            DM param_sens = optResult_lam_x(Slice(optResult_lam_x.size1() - varp.size1(), optResult_lam_x.size1()));

            for (auto &param : paramList) {
                param.foundValue = static_cast<double>(new_params(param.id));
                param.sensitivity = std::abs(static_cast<double>(param_sens(param.id)));
            }

            /* Time-reverse trajectory from optimization result, reshape, reverse time */
            DM est_traj_rev_vect = optResult_x(Slice(0, varx.size1()));
            DM est_traj_rev_woTime = DM::reshape(est_traj_rev_vect, dimx, DATA_POINTS);
            DM est_traj_woTime = flipDM(est_traj_rev_woTime, 2);

            /* Now: Rows: States, Columns: Time steps
             * Add original times as first row */
            DM est_traj_wTime = DM::vertcat({id_time, est_traj_woTime});

            /* Calculate fitting error (cost function) using time-reverse trajectory */
            double fitting_error_evaluated = static_cast<double>(fitting_error_func(DMVector{est_traj_rev_woTime})[0]);

            std::cout << "\n\nIteration " << iIt + 1 << " result:\n";

            /** Write estimated trajectory to file **/
            std::string est_traj_filepath = kite_params_out_dir + afterItStr + "_estimatedTrajectory.txt";
            write_trajectoryFile(est_traj_wTime, est_traj_filepath);

            /** Write new parameters to YAML file **/
            write_parameterFile(kite_params_in_filepath, paramList, kite_params_out_filepath);

            /** Visualize new parameters within their bounds and write to paramInfo file **/
            std::string paramInfo_filepath = kite_params_out_dir + afterItStr + "_paramInfo.txt";
            printWrite_parametersFound(paramList, fitting_error_evaluated, solve_status, paramInfo_filepath);

            /** Write optimization results to file for possible warm start **/
            write_optResultsFile(optResult_x, optResult_lam_x, optResult_lam_g, optResults0_filepath);


            /* Check on max iteration steps */
            if (i_optimization_failed >= n_optimizationFailed_max) {
                std::cout << "Optimization failed to many times. Stopping iterations on this sequence. \n";
                iIt = nIt;
                continue;
            }
        } /** End of ident iteration loop -------------------------------------------------------------------------- **/

    } /** End of sequence loop ------------------------------------------------------------------------------------- **/

}
