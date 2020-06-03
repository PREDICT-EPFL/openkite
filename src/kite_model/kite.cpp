#include "kite.h"

using namespace casadi;

namespace kite_utils {
KiteProperties LoadProperties(const std::string &filename) {
    //read YAML config file
    YAML::Node config = YAML::LoadFile(filename);

    //create properties object and fill in with data
    KiteProperties props;

    props.name = config["info"]["name"].as<std::string>();

    /** --------------------- **/
    /** Geometric parameters  **/
    /** --------------------- **/
    props.Geometry.wingSpan = config["geom"]["b"].as<double>();
    props.Geometry.mac = config["geom"]["c"].as<double>();
    props.Geometry.aspectRatio = config["geom"]["AR"].as<double>();
    props.Geometry.wingSurfaceArea = config["geom"]["S"].as<double>();
    //props.Geometry.TaperRatio = config["geometry"]["lam"].as<double>();
    //props.Geometry.HTailsurface = config["geometry"]["St"].as<double>();
    //props.Geometry.TailLeverArm = config["geometry"]["lt"].as<double>();
    //props.Geometry.FinSurfaceArea = config["geometry"]["Sf"].as<double>();
    //props.Geometry.FinLeverArm = config["geometry"]["lf"].as<double>();
    //props.Geometry.AerodynamicCenter = config["geometry"]["Xac"].as<double>();

    /** --------------------------- **/
    /** Mass and inertia parameters **/
    /** --------------------------- **/
    props.Inertia.mass = config["inertia"]["mass"].as<double>();
    props.Inertia.Ixx = config["inertia"]["Ixx"].as<double>();
    props.Inertia.Iyy = config["inertia"]["Iyy"].as<double>();
    props.Inertia.Izz = config["inertia"]["Izz"].as<double>();
    props.Inertia.Ixz = config["inertia"]["Ixz"].as<double>();

    /** ------------------------------- **/
    /** Aerodynamic parameters          **/
    /** ------------------------------- **/
    props.Aerodynamics.e_oswald = config["aero"]["e_oswald"].as<double>();

    props.Aerodynamics.CD0 = config["aero"]["CD0"].as<double>();


    /* AOA */
    props.Aerodynamics.CL0 = config["aero_aoa"]["CL0"].as<double>();
    props.Aerodynamics.CLa = config["aero_aoa"]["CLa"].as<double>();

    props.Aerodynamics.Cm0 = config["aero_aoa"]["Cm0"].as<double>();
    props.Aerodynamics.Cma = config["aero_aoa"]["Cma"].as<double>();

    /* Sideslip */
    props.Aerodynamics.CYb = config["aero_ss"]["CYb"].as<double>();

    props.Aerodynamics.Cl0 = config["aero_ss"]["Cl0"].as<double>();
    props.Aerodynamics.Clb = config["aero_ss"]["Clb"].as<double>();

    props.Aerodynamics.Cn0 = config["aero_ss"]["Cn0"].as<double>();
    props.Aerodynamics.Cnb = config["aero_ss"]["Cnb"].as<double>();


    /* Pitchrate */
    props.Aerodynamics.CLq = config["aero_rate_pitch"]["CLq"].as<double>();
    props.Aerodynamics.Cmq = config["aero_rate_pitch"]["Cmq"].as<double>();

    /* Rollrate */
    props.Aerodynamics.CYp = config["aero_rate_roll"]["CYp"].as<double>();
    props.Aerodynamics.Clp = config["aero_rate_roll"]["Clp"].as<double>();
    props.Aerodynamics.Cnp = config["aero_rate_roll"]["Cnp"].as<double>();

    /* Yawrate */
    props.Aerodynamics.CYr = config["aero_rate_yaw"]["CYr"].as<double>();
    props.Aerodynamics.Clr = config["aero_rate_yaw"]["Clr"].as<double>();
    props.Aerodynamics.Cnr = config["aero_rate_yaw"]["Cnr"].as<double>();


    /** ------------------------------ **/
    /** Aerodynamic effects of control **/
    /** ------------------------------ **/
    /* Elevator */
    props.Aerodynamics.CLde = config["aero_ctrl_elev"]["CLde"].as<double>();
    props.Aerodynamics.Cmde = config["aero_ctrl_elev"]["Cmde"].as<double>();

    /* Ailerons */
    props.Aerodynamics.Clda = config["aero_ctrl_ail"]["Clda"].as<double>();
    props.Aerodynamics.Cnda = config["aero_ctrl_ail"]["Cnda"].as<double>();

    /* Rudder */
    props.Aerodynamics.CYdr = config["aero_ctrl_rud"]["CYdr"].as<double>();
    props.Aerodynamics.Cldr = config["aero_ctrl_rud"]["Cldr"].as<double>();
    props.Aerodynamics.Cndr = config["aero_ctrl_rud"]["Cndr"].as<double>();


    /** ------------------------------ **/
    /** Actuator dynamics              **/
    /** ------------------------------ **/
    props.Actuators.TC_dER = config["actuators"]["TC_dER"].as<double>();
    props.Actuators.TC_dA = config["actuators"]["TC_dA"].as<double>();


    //props.Tether.length = config["tether"]["length"].as<double>();
    //props.Tether.Ks = config["tether"]["Ks"].as<double>();
    //props.Tether.Kd = config["tether"]["Kd"].as<double>();
    //props.Tether.rx = config["tether"]["rx"].as<double>();
    //props.Tether.ry = config["tether"]["ry"].as<double>();
    //props.Tether.rz = config["tether"]["rz"].as<double>();

    return props;
}

KiteProperties LoadMinimalProperties(const std::string &filename) {
    //read YAML config file
    YAML::Node config = YAML::LoadFile(filename);

    //create properties object and fill in with data
    KiteProperties props;

    props.name = config["info"]["name"].as<std::string>();

    /** --------------------- **/
    /** Geometric parameters  **/
    /** --------------------- **/
    props.Geometry.wingSpan = config["geom"]["b"].as<double>();
    props.Geometry.mac = config["geom"]["c"].as<double>();
    props.Geometry.aspectRatio = config["geom"]["AR"].as<double>();
    props.Geometry.wingSurfaceArea = config["geom"]["S"].as<double>();

    /** --------------------------- **/
    /** Mass and inertia parameters **/
    /** --------------------------- **/
    props.Inertia.mass = config["inertia"]["mass"].as<double>();
    props.Inertia.Ixx = config["inertia"]["Ixx"].as<double>();
    props.Inertia.Iyy = config["inertia"]["Iyy"].as<double>();
    props.Inertia.Izz = config["inertia"]["Izz"].as<double>();
    props.Inertia.Ixz = config["inertia"]["Ixz"].as<double>();

    /** ------------------------------- **/
    /** Aerodynamic parameters          **/
    /** ------------------------------- **/
    props.Aerodynamics.e_oswald = config["aero"]["e_oswald"].as<double>();

    props.Aerodynamics.CD0 = config["aero"]["CD0"].as<double>();


    /* AOA */
    props.Aerodynamics.CL0 = config["aero_aoa"]["CL0"].as<double>();
    props.Aerodynamics.CLa = config["aero_aoa"]["CLa"].as<double>();

    props.Aerodynamics.Cm0 = config["aero_aoa"]["Cm0"].as<double>();
    props.Aerodynamics.Cma = config["aero_aoa"]["Cma"].as<double>();

    /* Sideslip */
    props.Aerodynamics.CYb = config["aero_ss"]["CYb"].as<double>();

    props.Aerodynamics.Clb = config["aero_ss"]["Clb"].as<double>();

    props.Aerodynamics.Cnb = config["aero_ss"]["Cnb"].as<double>();

    /* Rollrate */
//        props.Aerodynamics.CYp = config["aero_rate_roll"]["CYp"].as<double>();
    props.Aerodynamics.Clp = config["aero_rate_roll"]["Clp"].as<double>();
    //   props.Aerodynamics.Cnp = config["aero_rate_roll"]["Cnp"].as<double>();

    /* Yawrate */
//        props.Aerodynamics.CYr = config["aero_rate_yaw"]["CYr"].as<double>();
    //    props.Aerodynamics.Clr = config["aero_rate_yaw"]["Clr"].as<double>();
    props.Aerodynamics.Cnr = config["aero_rate_yaw"]["Cnr"].as<double>();


    /** ------------------------------ **/
    /** Aerodynamic effects of control **/
    /** ------------------------------ **/
    /* Elevator */
//        props.Aerodynamics.CLde = config["aero_ctrl_elev"]["CLde"].as<double>();
    props.Aerodynamics.Cmde = config["aero_ctrl_elev"]["Cmde"].as<double>();

    /* Ailerons */
    props.Aerodynamics.Clda = config["aero_ctrl_ail"]["Clda"].as<double>();
//        props.Aerodynamics.Cnda = config["aero_ctrl_ail"]["Cnda"].as<double>();

    /* Rudder */
//        props.Aerodynamics.CYdr = config["aero_ctrl_rud"]["CYdr"].as<double>();
//        props.Aerodynamics.Cldr = config["aero_ctrl_rud"]["Cldr"].as<double>();
    props.Aerodynamics.Cndr = config["aero_ctrl_rud"]["Cndr"].as<double>();

    //props.Tether.length = config["tether"]["length"].as<double>();
    //props.Tether.Ks = config["tether"]["Ks"].as<double>();
    //props.Tether.Kd = config["tether"]["Kd"].as<double>();
    //props.Tether.rx = config["tether"]["rx"].as<double>();
    //props.Tether.ry = config["tether"]["ry"].as<double>();
    //props.Tether.rz = config["tether"]["rz"].as<double>();

    return props;
}


time_point get_time() {
    /** OS dependent */
#ifdef __APPLE__
    return std::chrono::system_clock::now();
#else
    return std::chrono::high_resolution_clock::now();
#endif
}

bool file_exists(const std::string &filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

DM read_from_file(const std::string &filename) {
    std::ifstream file(filename, std::ios::in);
    std::vector<double> vec;
    if (!file.fail()) {
        double x;
        while (file >> x) {
            vec.push_back(x);
        }
        return DM({vec});
    } else {
        std::cout << "Could not open : " << filename << " data file \n";
        file.clear();
        return DM({});
    }
}

void write_to_file(const std::string &filename, const DM &data) {
    std::ofstream data_file(filename, std::ios::out);
    std::vector<double> vec = data.nonzeros();

    /** solution */
    if (!data_file.fail()) {
        for (std::vector<double>::iterator it = vec.begin(); it != vec.end(); ++it) {
            data_file << (*it) << " ";
        }
        data_file << "\n";
    }
    data_file.close();
}
}


/* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder */
template<typename W, typename GEN, typename DLO, typename DLA, typename AIL, typename ELV, typename RUD>
void KiteDynamics::getModel(GEN &g, GEN &rho,
                            W &windFrom, W &windSpeed,
                            GEN &b, GEN &c, GEN &AR, GEN &S,
                            GEN &Mass, GEN &Ixx, GEN &Iyy, GEN &Izz, GEN &Ixz,

                            DLO &e_oswald,
                            DLO &CD0,


                            DLO &CL0,
                            DLO &CLa,

                            DLO &Cm0,
                            DLO &Cma,


                            DLA &CYb,

                            GEN &Cl0,
                            DLA &Clb,

                            GEN &Cn0,
                            DLA &Cnb,


                            DLO &CLq,
                            DLO &Cmq,

                            DLA &CYp,
                            DLA &Clp,
                            DLA &Cnp,

                            DLA &CYr,
                            DLA &Clr,
                            DLA &Cnr,


                            ELV &CLde,
                            ELV &Cmde,

                            AIL &Clda,
                            AIL &Cnda,

                            RUD &CYdr,
                            RUD &Cldr,
                            RUD &Cndr,


                            ELV &TC_dE,
                            RUD &TC_dR,
                            AIL &TC_dA,

                            SX &v, SX &w, SX &r, SX &q,
                            SX &F_thr0, SX &dE, SX &dR, SX &dA,

                            SX &F_thr0_cmd, SX &dE_cmd, SX &dR_cmd, SX &dA_cmd,

                            SX &v_dot, SX &w_dot, SX &r_dot, SX &q_dot,
                            SX &F_thr0_dot, SX &dE_dot, SX &dR_dot, SX &dA_dot,

                            SX &Va_pitot, SX &Va, SX &alpha, SX &beta,
                            SX &b_F_aero, SX &b_F_thrust, bool teth_ON, SX &b_F_tether) {
    /** Start of model ============================================================================================= **/
    /* Aircraft Inertia Matrix */
    auto J = SX::diag(SX::vertcat({Ixx, Iyy, Izz}));
    J(0, 2) = Ixz;
    J(2, 0) = Ixz;

    /** State variables **/
    v = SX::sym("v", 3); // Linear velocity of the CoG (body frame) [m/s]
    w = SX::sym("w", 3); // Angular velocity (body frame) [rad/s]
    r = SX::sym("r", 3); // Position of the CoG (geodetic (NED) frame) [m]
    q = SX::sym("q", 4); // Rotation from geodetic (NED) to body frame. Transforms a body vector to ned. = q_gb
    SX q_bg = kmath::quat_inverse(q);

    /** Control commands **/
    F_thr0_cmd = SX::sym("F_thr0_cmd");
    dE_cmd = SX::sym("dE_cmd");
    dR_cmd = SX::sym("dR_cmd");
    dA_cmd = SX::sym("dA_cmd");

    /** Control positions (states) **/
    F_thr0 = SX::sym("F_thr0");   // Static propeller thrust along body frame x axis [N]
    dE = SX::sym("dE"); // Elevator deflection [positive causing negative pitch movement (nose down)] [rad]
    dR = SX::sym("dR"); // Rudder deflection [positive causing negative yaw movement (nose left)] [rad]
    dA = SX::sym("dA"); // Aileron deflection [positive causing negative roll movement (right wing up)] [rad]

    /** Aerodynamic (Wind) frame **/
    /* Wind velocity in body frame */
    SX g_vW = windSpeed * SX::vertcat({-cos(windFrom), -sin(windFrom), 0});
    SX b_vW = kmath::quat_transform(q_bg, g_vW);

    /* Apparent velocity in body frame, airspeed */
    SX b_va = v - b_vW;

    /* Aerodynamic variables (Airspeed, angle of attack, side slip angle) */
    Va = SX::norm_2(b_va);
    alpha = atan(b_va(2) / (b_va(0))); // + 1e-4));
    beta = asin(b_va(1) / (Va)); // + 1e-4));

    /* Measured airspeed component (pitot tube orientation dependent) */
    SX r_sens = SX::vertcat({0.11, 0.22, -0.05});
    SX b_va_meas = b_va + SX::cross(r_sens, w);
    /* At fast body yawrates, the body rotation (thus pitot tube is faster than the CoG) is clearly negligible
     * under the effect of the pitot measurement direction being rotated out of the airflow (sideslip) */

    /* Pitot tube is oriented about 5 degress above body x axes */
    SX q_sens_b = kmath::T2quat(5.0 * M_PI / 180.0);
    SX sens_va = kmath::quat_transform(q_sens_b, b_va_meas);
    Va_pitot = sens_va(0);

    /** ---------------------------------------------------------- **/
    /** Aerodynamic Forces and Moments in aerodynamic (wind) frame **/
    /** ---------------------------------------------------------- **/
//    SX q_ba = kmath::quat_multiply(kmath::T2quat(alpha), kmath::T2quat(-beta));
//    SX w_ = w;

    /* XFLR5 gives coefficients in stability frame
     * To get from stability axis to body, rotate by aoa */
    SX q_bs = kmath::T2quat(alpha);
    /* Body rates in stability frame (pitch rate is equal in both frames) */
    SX q_sb = kmath::quat_inverse(q_bs);
    SX w_ = kmath::quat_transform(q_sb, w);

    SX dyn_press = 0.5 * rho * Va * Va;
    SX CL = CL0 + CLa * alpha + CLq * c / (2.0 * Va) * w_(1) + CLde * dE;
    SX CD = (CD0 + CL * CL / (pi * e_oswald * AR));

    /** Forces in x, y, z directions: -Drag, Side force, -Lift **/
    SX LIFT = dyn_press * S * CL;
    SX DRAG = dyn_press * S * CD;
    SX SF = dyn_press * S * (CYb * beta +
                             b / (2.0 * Va) * (CYp * w_(0) + CYr * w_(2)) +
                             CYdr * dR);

    SX Faero = SX::vertcat({-DRAG, SF, -LIFT});

    /** Moments about x, y, z axes: L, M, N **/
    SX L = dyn_press * S * b * (Cl0 + Clb * beta +
                                b / (2.0 * Va) * (Clp * w_(0) + Clr * w_(2)) +
                                Clda * dA + Cldr * dR);

    SX M = dyn_press * S * c * (Cm0 + Cma * alpha +
                                c / (2.0 * Va) * Cmq * w_(1) +
                                Cmde * dE);

    SX N = dyn_press * S * b * (Cn0 + Cnb * beta +
                                b / (2.0 * Va) * (Cnp * w_(0) + Cnr * w_(2)) +
                                Cnda * dA + Cndr * dR);

    SX Maero = SX::vertcat({L, M, N});

    /** Aerodynamic Forces and Moments in body frame **/
//    b_F_aero = kmath::quat_transform(q_ba, Faero);
    b_F_aero = kmath::quat_transform(q_bs, Faero);

//    auto b_Maero = kmath::quat_transform(q_ba, Maero);
    auto b_Maero = kmath::quat_transform(q_bs, Maero);
//    auto b_Maero =  Maero;

    /** ---------------------------------------- **/
    /** Gravitation, Thrust, Tether (body frame) **/
    /** ---------------------------------------- **/
    /** Gravitational acceleration **/
    SX b_g = kmath::quat_transform(q_bg, SX::vertcat({0, 0, g}));

    /** Propeller thrust **/
    const double p1 = -0.014700129;
    const double p2 = -0.00083119;
    b_F_thrust = F_thr0 * (p2 * Va * Va + p1 * Va + 1.0) * SX::vertcat({1, 0, 0});
    // F_thr0 is the static thrust (at zero airspeed)

    /** Tether force and moment **/
    if (teth_ON) {
        SX dist = SX::norm_2(r);
        const double tethLen = 120.0;
        const double tethDiameter = 0.3e-3;
        const double tethCrossArea = M_PI_4 * tethDiameter * tethDiameter;
        const double tethDensity = 1.15;
        const double tethMassPerMeter = tethCrossArea * 1.0 * tethDensity;
        const double tethE = 9.05e9;

        const double c_orth = 1.2;

        /* Weight_tether */
        //SX decl = SX::atan(SX::norm_2(r(Slice(0, 2))) / -r(2));
        //SX g_Wteth = 0.5 * (1 + cos(decl)) * tethLen * tethMassPerMeter * SX::vertcat({0, 0, g}); // probably invalid
        SX g_Wteth = tethLen * tethMassPerMeter * SX::vertcat({0, 0, g});
        SX b_Wteth = kmath::quat_transform(q_bg, g_Wteth);

        /* Drag_tether */
        SX q_ba = kmath::quat_multiply(kmath::T2quat(alpha), kmath::T3quat(-beta));

        SX lat = SX::atan2(-r(0), SX::norm_2(r(Slice(1, 3))));
        SX lon = SX::atan2(-r(1), -r(2));
        SX q_lg = kmath::quat_multiply(kmath::T2quat(-lat), kmath::T1quat(lon));
        SX q_gl = kmath::quat_inverse(q_lg);
        SX q_lb = kmath::quat_multiply(q_lg, q);

        SX l_va = kmath::quat_transform(q_lb, b_va);
        SX l_va_proj = SX::vertcat({l_va(0), l_va(1), 0});

        SX q_bl = kmath::quat_multiply(q_bg, q_gl);
        SX f_va_proj = kmath::quat_transform(q_bl, l_va_proj);

        SX b_Dteth = 0.125 * tethDensity * -f_va_proj * SX::norm_2(f_va_proj) * c_orth * tethDiameter * tethLen;

        /* LonForce_tether */
        SX tethElongation = (dist - tethLen) / tethLen;
        SX g_lonFteth = -r / dist * SX::fmax(0, (tethE * tethCrossArea) * (tethElongation + 0.005));
        SX b_lonFteth = kmath::quat_transform(q_bg, g_lonFteth);
        b_F_tether = b_lonFteth + b_Dteth + b_Wteth;
    } else {
        b_F_tether = SX::vertcat({0, 0, 0});
    }

    SX tether_outlet_position = SX::vertcat({0, 0, 0.05});
    SX b_Mtether = SX::cross(tether_outlet_position, b_F_tether);

    /** ----------------------------- **/
    /** Motion equations (body frame) **/
    /** ----------------------------- **/
    /** Linear motion equation **/
    v_dot = (b_F_aero + b_F_thrust + b_F_tether) / Mass + b_g - SX::cross(w, v);

    /** Angular motion equation **/
    w_dot = SX::mtimes(SX::inv(J), (b_Maero + b_Mtether - SX::cross(w, SX::mtimes(J, w))));

    /** ------------------------------------ **/
    /** Kinematic Equations (geodetic frame) **/
    /** ------------------------------------ **/
    /** Translation: Aircraft position derivative **/
    r_dot = kmath::quat_transform(q, v);                            // q = q_gb, v (body frame)

    /** Rotation: Aircraft attitude derivative **/
    double lambda = -5;
    q_dot = 0.5 * kmath::quat_multiply(q, SX::vertcat({0, w}))      // q = q_gb, w = omega (body frame)
            + 0.5 * lambda * q * (SX::dot(q, q) - 1);               // Quaternion norm stabilization term,
    // as in Gros: 'Baumgarte Stabilisation over the SO(3) Rotation Group for Control',
    // improved: lambda negative and SX::dot(q, q) instead of lambda positive and 1/SX::dot(q, q).

    /** ------------------------------------ **/
    /** Actuator dynamics **/
    /** ------------------------------------ **/
    const double TC_thr = 1.0 / 6.3739;               // Full thrust builds up in 1 second
    // const double TC_dE_ = 0.11 / (60.0 * M_PI / 180.0);
    // const double TC_dR_ = 0.11 / (60.0 * M_PI / 180.0);
    // const double TC_dA_ = 0.12 / (40.0 * M_PI / 180.0);

    F_thr0_dot = (F_thr0_cmd - F_thr0) / TC_thr;
    dE_dot = (dE_cmd - dE) / TC_dE;
    dR_dot = (dR_cmd - dR) / TC_dR;
    dA_dot = (dA_cmd - dA) / TC_dA;

    /** End of dynamics model ====================================================================================== **/
}

KiteDynamics::KiteDynamics(const KiteProperties &kiteProps, const AlgorithmProperties &AlgoProps,
                           const bool teth_ON) {

    /** enviromental constants */
    double g = 9.806;
    double rho = kiteProps.atmosphere.airDensity;

    /** --------------------- **/
    /** Geometric parameters  **/
    /** --------------------- **/
    double b = kiteProps.Geometry.wingSpan;
    double c = kiteProps.Geometry.mac;
    double AR = kiteProps.Geometry.aspectRatio;
    double S = kiteProps.Geometry.wingSurfaceArea;

    /** --------------------------- **/
    /** Mass and inertia parameters **/
    /** --------------------------- **/
    double Mass = kiteProps.Inertia.mass;
    double Ixx = kiteProps.Inertia.Ixx;
    double Iyy = kiteProps.Inertia.Iyy;
    double Izz = kiteProps.Inertia.Izz;
    double Ixz = kiteProps.Inertia.Ixz;

    /** ------------------------------- **/
    /** Aerodynamic parameters          **/
    /** ------------------------------- **/
    double e_oswald = kiteProps.Aerodynamics.e_oswald;
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

    /** ------------------------------ **/
    /** Actuator dynamics              **/
    /** ------------------------------ **/
    double TC_dE = kiteProps.Actuators.TC_dER;
    double TC_dR = kiteProps.Actuators.TC_dER;
    double TC_dA = kiteProps.Actuators.TC_dA;

    /** ------------------------------ **/
    /**        Tether parameters       **/
    /** ------------------------------ **/
    //    double Lt = KiteProps.Tether.length;
    //    double Ks = KiteProps.Tether.Ks;
    //    double Kd = KiteProps.Tether.Kd;
    //    double rx = KiteProps.Tether.rx;
    //    double ry = KiteProps.Tether.ry;
    //    double rz = KiteProps.Tether.rz;


    SX v, w, r, q, F_thr0, dE, dR, dA;
    SX F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd;
    SX v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot;
    SX Va_pitot, Va, alpha, beta;
    SX b_F_aero, b_F_thrust, b_F_tether;

    double windFrom = kiteProps.atmosphere.WindFrom;
    double windSpeed = kiteProps.atmosphere.WindSpeed;

    /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
     *        W,     GEN,   DLO,     DLA,   AIL,    ELV,    RUD*/
    getModel<double, double, double, double, double, double, double>(
            g, rho,
            windFrom, windSpeed,
            b, c, AR, S,
            Mass, Ixx, Iyy, Izz, Ixz,

            e_oswald, CD0,

            CL0, CLa,
            Cm0, Cma,

            CYb,
            Cl0, Clb,
            Cn0, Cnb,

            CLq, Cmq,
            CYp, Clp, Cnp,
            CYr, Clr, Cnr,

            CLde, Cmde,
            Clda, Cnda,
            CYdr, Cldr, Cndr,

            TC_dE, TC_dR, TC_dA,

            v, w, r, q, F_thr0, dE, dR, dA,
            F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd,
            v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot,

            Va_pitot, Va, alpha, beta,
            b_F_aero, b_F_thrust, teth_ON, b_F_tether);

    SX control_cmd = SX::vertcat({F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd});

    /** Complete dynamics of the Kite */
    auto state = SX::vertcat({v, w, r, q, F_thr0, dE, dR, dA});
    auto dynamics = SX::vertcat({v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot});

    SX aero_values = SX::vertcat({Va, alpha, beta});

    Function dyn_func = Function("dynamics", {state, control_cmd}, {dynamics});
    Function airspeedMeas_func = Function("airspeed", {state}, {Va_pitot});
    Function aeroValues_func = Function("aero_out", {state}, {aero_values});
    Function specNongravForce_func = Function("spec_nongrav_force", {state, control_cmd},
                                              {(b_F_aero + b_F_thrust + b_F_tether) / Mass});
    Function specTethForce_func = Function("specTethForce", {state, control_cmd}, {b_F_tether / Mass});

    /** compute dynamics state Jacobian */
    SX d_jacobian = SX::jacobian(dynamics, state);
    Function dyn_jac = Function("dyn_jacobian", {state, control_cmd}, {d_jacobian});

    AeroDynamics = Function("Aero", {state, control_cmd}, {b_F_aero});
    /** define RK4 integrator scheme */
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control_cmd.size1());
    SX dT = SX::sym("dT");

    /** get symbolic expression for RK4 integrator */
    SX sym_integrator = kmath::rk4_symbolic(X, U, dyn_func, dT);
    Function RK4_INT = Function("RK4", {X, U, dT}, {sym_integrator});

    /** CVODES integrator */
    /** @todo: make smarter initialisation of integrator */
    double h = AlgoProps.sampling_time;
    SXDict ode = {{"x",   state},
                  {"p",   control_cmd},
                  {"ode", dynamics}};
    Dict opts = {{"tf", h}};
    Function CVODES_INT = integrator("CVODES_INT", "cvodes", ode, opts);

    /** assign class atributes */
    this->State = state;
    this->Control = control_cmd;
    this->SymDynamics = dynamics;
    this->SymIntegartor = sym_integrator;
    this->SymJacobian = d_jacobian;

    this->NumDynamics = dyn_func;
    this->NumAirspeedMeasured = airspeedMeas_func;
    this->NumAeroValues = aeroValues_func;
    this->NumSpecNongravForce = specNongravForce_func;
    this->NumSpecTethForce = specTethForce_func;
    this->NumJacobian = dyn_jac;

    /** return integrator function */
    if (AlgoProps.Integrator == IntType::CVODES)
        this->NumIntegrator = CVODES_INT;
    else
        this->NumIntegrator = RK4_INT;

}

KiteDynamics::KiteDynamics(const KiteProperties &kiteProps, const AlgorithmProperties &AlgoProps,
                           const kite_utils::IdentMode &identMode) {

    /** enviromental constants */
    double g = 9.806;
    double rho = kiteProps.atmosphere.airDensity;
    double windFrom = kiteProps.atmosphere.WindFrom;
    double windSpeed = kiteProps.atmosphere.WindSpeed;

    /** --------------------- **/
    /** Geometric parameters  **/
    /** --------------------- **/
    double b = kiteProps.Geometry.wingSpan;
    double c = kiteProps.Geometry.mac;
    double AR = kiteProps.Geometry.aspectRatio;
    double S = kiteProps.Geometry.wingSurfaceArea;

    /** --------------------------- **/
    /** Mass and inertia parameters **/
    /** --------------------------- **/
    double Mass = kiteProps.Inertia.mass;
    double Ixx = kiteProps.Inertia.Ixx;
    double Iyy = kiteProps.Inertia.Iyy;
    double Izz = kiteProps.Inertia.Izz;
    double Ixz = kiteProps.Inertia.Ixz;

    SX v, w, r, q, F_thr0, dE, dR, dA;
    SX F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd;
    SX v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot;
    SX Va_pitot, Va, alpha, beta;
    SX b_F_aero, b_F_thrust, b_F_tether;

    SX params;
    SX control_cmd;

    if (identMode == kite_utils::IdentMode::PITCH) {
        /** LONGITUDINAL IDENTIFICATION PARAMETERS ------------------------------------------------------------------ */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
        //double e_oswald = kiteProps.Aerodynamics.e_oswald;
        SX e_oswald = SX::sym("e_oswald");
        SX CD0 = SX::sym("CD0");

        /* AOA */
        SX CL0 = SX::sym("CL0");
        SX CLa = SX::sym("CLa");

        SX Cm0 = SX::sym("Cm0");
        SX Cma = SX::sym("Cma");

        /* Sideslip */
        double CYb = kiteProps.Aerodynamics.CYb;

        double Cl0 = kiteProps.Aerodynamics.Cl0;
        double Clb = kiteProps.Aerodynamics.Clb;

        double Cn0 = kiteProps.Aerodynamics.Cn0;
        double Cnb = kiteProps.Aerodynamics.Cnb;


        /* Pitchrate */
        SX CLq = SX::sym("CLq");
        SX Cmq = SX::sym("Cmq");

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
        SX CLde = SX::sym("CLde");
        SX Cmde = SX::sym("Cmde");

        /* Ailerons */
        double Clda = kiteProps.Aerodynamics.Clda;
        double Cnda = kiteProps.Aerodynamics.Cnda;

        /* Rudder */
        double CYdr = kiteProps.Aerodynamics.CYdr;
        double Cldr = kiteProps.Aerodynamics.Cldr;
        double Cndr = kiteProps.Aerodynamics.Cndr;


        /** ------------------------------ **/
        /** Actuator dynamics              **/
        /** ------------------------------ **/
        SX TC_dE = SX::sym("TC_dE");
        double TC_dR = kiteProps.Actuators.TC_dER;
        double TC_dA = kiteProps.Actuators.TC_dA;

        params = SX::vertcat({params,
                              e_oswald,
                              CD0,

                              CL0,
                              CLa,

                              Cm0,
                              Cma,

                              CLq,
                              Cmq,

                              CLde,
                              Cmde,

                              TC_dE
                             }); // 11 longitudinal parameters

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W, GEN,   DLO, DLA,   AIL,    ELV, RUD*/
        getModel<double, double, SX, double, double, SX, double>(
                g, rho,
                windFrom, windSpeed,
                b, c, AR, S,
                Mass, Ixx, Iyy, Izz, Ixz,

                e_oswald, CD0,

                CL0, CLa,
                Cm0, Cma,

                CYb,
                Cl0, Clb,
                Cn0, Cnb,

                CLq, Cmq,
                CYp, Clp, Cnp,
                CYr, Clr, Cnr,

                CLde, Cmde,
                Clda, Cnda,
                CYdr, Cldr, Cndr,

                TC_dE, TC_dR, TC_dA,

                v, w, r, q, F_thr0, dE, dR, dA,
                F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd,
                v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot,

                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, false, b_F_tether);

        control_cmd = SX::vertcat({F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd});


    } else if (identMode == kite_utils::IdentMode::ROLL) {
        /** LATERAL IDENTIFICATION PARAMETERS ----------------------------------------------------------------------- */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
        double e_oswald = kiteProps.Aerodynamics.e_oswald;
        double CD0 = kiteProps.Aerodynamics.CD0;

        /* AOA */
        double CL0 = kiteProps.Aerodynamics.CL0;
        double CLa = kiteProps.Aerodynamics.CLa;

        double Cm0 = kiteProps.Aerodynamics.Cm0;
        double Cma = kiteProps.Aerodynamics.Cma;

        /* Sideslip */
        SX CYb = SX::sym("CYb");

        double Cl0 = kiteProps.Aerodynamics.Cl0;
        SX Clb = SX::sym("Clb");

        double Cn0 = kiteProps.Aerodynamics.Cn0;
        SX Cnb = SX::sym("Cnb");


        /* Pitchrate */
        double CLq = kiteProps.Aerodynamics.CLq;
        double Cmq = kiteProps.Aerodynamics.Cmq;

        /* Rollrate */
        SX CYp = SX::sym("CYp");
        SX Clp = SX::sym("Clp");
        SX Cnp = SX::sym("Cnp");

        /* Yawrate */
        SX CYr = SX::sym("CYr");
        SX Clr = SX::sym("Clr");
        SX Cnr = SX::sym("Cnr");


        /** ------------------------------ **/
        /** Aerodynamic effects of control **/
        /** ------------------------------ **/
        /* Elevator */
        double CLde = kiteProps.Aerodynamics.CLde;
        double Cmde = kiteProps.Aerodynamics.Cmde;

        /* Ailerons */
        SX Clda = SX::sym("Clda");
        SX Cnda = SX::sym("Cnda");

        /* Rudder */
        double CYdr = kiteProps.Aerodynamics.CYdr;
        double Cldr = kiteProps.Aerodynamics.Cldr;
        double Cndr = kiteProps.Aerodynamics.Cndr;


        /** ------------------------------ **/
        /** Actuator dynamics              **/
        /** ------------------------------ **/
        double TC_dE = kiteProps.Actuators.TC_dER;
        double TC_dR = kiteProps.Actuators.TC_dER;
        SX TC_dA = SX::sym("TC_dA");

        params = SX::vertcat({params,
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
                              Cnda,

                              TC_dA
                             }); // 13 lateral parameters (Aileron control only)

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W,     GEN,    DLO,   DLA, AIL, ELV,   RUD*/
        getModel<double, double, double, SX, SX, double, double>(
                g, rho,
                windFrom, windSpeed,
                b, c, AR, S,
                Mass, Ixx, Iyy, Izz, Ixz,

                e_oswald, CD0,

                CL0, CLa,
                Cm0, Cma,

                CYb,
                Cl0, Clb,
                Cn0, Cnb,

                CLq, Cmq,
                CYp, Clp, Cnp,
                CYr, Clr, Cnr,

                CLde, Cmde,
                Clda, Cnda,
                CYdr, Cldr, Cndr,

                TC_dE, TC_dR, TC_dA,

                v, w, r, q, F_thr0, dE, dR, dA,
                F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd,
                v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot,
                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, false, b_F_tether);

        control_cmd = SX::vertcat({F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd});

    } else if (identMode == kite_utils::IdentMode::YAW) {
        /** LATERAL IDENTIFICATION PARAMETERS ----------------------------------------------------------------------- */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
        double e_oswald = kiteProps.Aerodynamics.e_oswald;
        double CD0 = kiteProps.Aerodynamics.CD0;

        /* AOA */
        double CL0 = kiteProps.Aerodynamics.CL0;
        double CLa = kiteProps.Aerodynamics.CLa;

        double Cm0 = kiteProps.Aerodynamics.Cm0;
        double Cma = kiteProps.Aerodynamics.Cma;

        /* Sideslip */
        SX CYb = SX::sym("CYb");

        double Cl0 = kiteProps.Aerodynamics.Cl0;
        SX Clb = SX::sym("Clb");

        double Cn0 = kiteProps.Aerodynamics.Cn0;
        SX Cnb = SX::sym("Cnb");


        /* Pitchrate */
        double CLq = kiteProps.Aerodynamics.CLq;
        double Cmq = kiteProps.Aerodynamics.Cmq;

        /* Rollrate */
        SX CYp = SX::sym("CYp");
        SX Clp = SX::sym("Clp");
        SX Cnp = SX::sym("Cnp");

        /* Yawrate */
        SX CYr = SX::sym("CYr");
        SX Clr = SX::sym("Clr");
        SX Cnr = SX::sym("Cnr");


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
        SX CYdr = SX::sym("CYdr");
        SX Cldr = SX::sym("Cldr");
        SX Cndr = SX::sym("Cndr");


        /** ------------------------------ **/
        /** Actuator dynamics              **/
        /** ------------------------------ **/
        double TC_dE = kiteProps.Actuators.TC_dER;
        SX TC_dR = SX::sym("TC_dR");
        double TC_dA = kiteProps.Actuators.TC_dA;

        params = SX::vertcat({params,
                              CYb,

                              Clb,

                              Cnb,

                              CYp,
                              Clp,
                              Cnp,

                              CYr,
                              Clr,
                              Cnr,

                              CYdr,
                              Cldr,
                              Cndr,

                              TC_dR
                             }); // 14 lateral parameters (Rudder control only)

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W,     GEN,    DLO,   DLA, AIL,    ELV,   RUD*/
        getModel<double, double, double, SX, double, double, SX>(
                g, rho,
                windFrom, windSpeed,
                b, c, AR, S,
                Mass, Ixx, Iyy, Izz, Ixz,

                e_oswald, CD0,

                CL0, CLa,
                Cm0, Cma,

                CYb,
                Cl0, Clb,
                Cn0, Cnb,

                CLq, Cmq,
                CYp, Clp, Cnp,
                CYr, Clr, Cnr,

                CLde, Cmde,
                Clda, Cnda,
                CYdr, Cldr, Cndr,

                TC_dE, TC_dR, TC_dA,

                v, w, r, q, F_thr0, dE, dR, dA,
                F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd,
                v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot,

                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, false, b_F_tether);

        control_cmd = SX::vertcat({F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd});

    } else if (identMode == kite_utils::IdentMode::COMPLETE) {
        /** COMPLETE IDENTIFICATION PARAMETERS ------------------------------------------------------------------ */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
//        double e_oswald = kiteProps.Aerodynamics.e_oswald;
        SX e_oswald = SX::sym("e_oswald");
        SX CD0 = SX::sym("CD0");

        /* AOA */
        SX CL0 = SX::sym("CL0");
        SX CLa = SX::sym("CLa");

        SX Cm0 = SX::sym("Cm0");
        SX Cma = SX::sym("Cma");

        /* Sideslip */
        SX CYb = SX::sym("CYb");

        double Cl0 = kiteProps.Aerodynamics.Cl0;
        SX Clb = SX::sym("Clb");

        double Cn0 = kiteProps.Aerodynamics.Cn0;
        SX Cnb = SX::sym("Cnb");


        /* Pitchrate */
        SX CLq = SX::sym("CLq");
        SX Cmq = SX::sym("Cmq");

        /* Rollrate */
        SX CYp = SX::sym("CYp");
        SX Clp = SX::sym("Clp");
        SX Cnp = SX::sym("Cnp");

        /* Yawrate */
        SX CYr = SX::sym("CYr");
        SX Clr = SX::sym("Clr");
        SX Cnr = SX::sym("Cnr");


        /** ------------------------------ **/
        /** Aerodynamic effects of control **/
        /** ------------------------------ **/
        /* Elevator */
        SX CLde = SX::sym("CLde");
        SX Cmde = SX::sym("Cmde");

        /* Ailerons */
        SX Clda = SX::sym("Clda");
        SX Cnda = SX::sym("Cnda");

        /* Rudder */
        SX CYdr = SX::sym("CYdr");
        SX Cldr = SX::sym("Cldr");
        SX Cndr = SX::sym("Cndr");

        /** ------------------------------ **/
        /** Actuator dynamics              **/
        /** ------------------------------ **/
        SX TC_dE = SX::sym("TC_dE");
        SX TC_dR = SX::sym("TC_dR");
        SX TC_dA = SX::sym("TC_dA");

        params = SX::vertcat({params,
                              e_oswald,
                              CD0,

                              CL0,
                              CLa,

                              Cm0,
                              Cma,


                              CYb,

                              Clb,

                              Cn0,
                              Cnb,


                              CLq,
                              Cmq,

                              CYp,
                              Clp,
                              Cnp,

                              CYr,
                              Clr,
                              Cnr,


                              CLde,
                              Cmde,

                              Clda,
                              Cnda,

                              CYdr,
                              Cldr,
                              Cndr,

                              TC_dE,
                              TC_dR,
                              TC_dA
                             }); // 27 parameters

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W,     GEN,   DLO,     DLA,   AIL,    ELV,    RUD*/
        getModel<double, double, SX, SX, SX, SX, SX>(
                g, rho,
                windFrom, windSpeed,
                b, c, AR, S,
                Mass, Ixx, Iyy, Izz, Ixz,

                e_oswald, CD0,

                CL0, CLa,
                Cm0, Cma,

                CYb,
                Cl0, Clb,
                Cn0, Cnb,

                CLq, Cmq,
                CYp, Clp, Cnp,
                CYr, Clr, Cnr,

                CLde, Cmde,
                Clda, Cnda,
                CYdr, Cldr, Cndr,

                TC_dE, TC_dR, TC_dA,

                v, w, r, q, F_thr0, dE, dR, dA,
                F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd,
                v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot,

                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, false, b_F_tether);

        control_cmd = SX::vertcat({F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd});

    } else if (identMode == kite_utils::IdentMode::YR or identMode == kite_utils::IdentMode::YRb) {
        /** Yaw Roll IDENTIFICATION PARAMETERS ------------------------------------------------------------------ */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
        double e_oswald = kiteProps.Aerodynamics.e_oswald;
        double CD0 = kiteProps.Aerodynamics.CD0;

        /* AOA */
        double CL0 = kiteProps.Aerodynamics.CL0;
        double CLa = kiteProps.Aerodynamics.CLa;

        double Cm0 = kiteProps.Aerodynamics.Cm0;
        double Cma = kiteProps.Aerodynamics.Cma;

        /* Sideslip */
        SX CYb = SX::sym("CYb");

        double Cl0 = kiteProps.Aerodynamics.Cl0;
        SX Clb = SX::sym("Clb");

        double Cn0 = kiteProps.Aerodynamics.Cn0;
        SX Cnb = SX::sym("Cnb");

        /* Pitchrate */
        double CLq = kiteProps.Aerodynamics.CLq;
        double Cmq = kiteProps.Aerodynamics.Cmq;

        /* Rollrate */
        SX CYp = SX::sym("CYp");
        SX Clp = SX::sym("Clp");
        SX Cnp = SX::sym("Cnp");

        /* Yawrate */
        SX CYr = SX::sym("CYr");
        SX Clr = SX::sym("Clr");
        SX Cnr = SX::sym("Cnr");


        /** ------------------------------ **/
        /** Aerodynamic effects of control **/
        /** ------------------------------ **/
        /* Elevator */
        double CLde = kiteProps.Aerodynamics.CLde;
        double Cmde = kiteProps.Aerodynamics.Cmde;

        /* Ailerons */
        SX Clda = SX::sym("Clda");
        SX Cnda = SX::sym("Cnda");

        /* Rudder */
        SX CYdr = SX::sym("CYdr");
        SX Cldr = SX::sym("Cldr");
        SX Cndr = SX::sym("Cndr");


        /** ------------------------------ **/
        /** Actuator dynamics              **/
        /** ------------------------------ **/
        double TC_dE = kiteProps.Actuators.TC_dER;
        SX TC_dR = SX::sym("TC_dR");
        SX TC_dA = SX::sym("TC_dA");

        params = SX::vertcat({params,
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
                              Cnda,

                              CYdr,
                              Cldr,
                              Cndr,

                              TC_dE,
                              TC_dR,
                              TC_dA
                             }); // 16 parameters

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W,     GEN,   DLO,    DLA, AIL,  ELV,  RUD*/
        getModel<double, double, double, SX, SX, double, SX>(
                g, rho,
                windFrom, windSpeed,
                b, c, AR, S,
                Mass, Ixx, Iyy, Izz, Ixz,

                e_oswald, CD0,

                CL0, CLa,
                Cm0, Cma,

                CYb,
                Cl0, Clb,
                Cn0, Cnb,

                CLq, Cmq,
                CYp, Clp, Cnp,
                CYr, Clr, Cnr,

                CLde, Cmde,
                Clda, Cnda,
                CYdr, Cldr, Cndr,

                TC_dE, TC_dR, TC_dA,

                v, w, r, q, F_thr0, dE, dR, dA,
                F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd,
                v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot,

                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, false, b_F_tether);

        control_cmd = SX::vertcat({F_thr0_cmd, dE_cmd, dR_cmd, dA_cmd});
    }

//    /** ------------------------------ **/
//    /**        Tether parameters       **/
//    /** ------------------------------ **/
//    double Ks = kiteProps.Tether.Ks;
//    double Kd = kiteProps.Tether.Kd;
//    double Lt = kiteProps.Tether.length;
//    double rx = kiteProps.Tether.rx;
// //   double ry = kiteProps.Tether.ry;
//    double rz = kiteProps.Tether.rz;
//
//    params = SX::vertcat({params,
//                                 //                               Ks, Kd, Lt, rx, rz,
//                         }); // 5 Tether parameters

    /** Complete dynamics of the Kite */
    auto state = SX::vertcat({v, w, r, q, F_thr0, dE, dR, dA});
    auto dynamics = SX::vertcat({v_dot, w_dot, r_dot, q_dot, F_thr0_dot, dE_dot, dR_dot, dA_dot});

    Function dyn_func = Function("dynamics", {state, control_cmd, params}, {dynamics});
    Function specNongravForce_func = Function("spec_nongrav_force", {state, control_cmd},
                                              {(b_F_aero + b_F_thrust) / Mass});

    /** compute dynamics state Jacobian */
    SX d_jacobian = SX::jacobian(dynamics, state);
    Function dyn_jac = Function("dyn_jacobian", {state, control_cmd, params}, {d_jacobian});

    /** define RK4 integrator scheme */
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control_cmd.size1());
    SX dT = SX::sym("dT");

    /** get symbolic expression for RK4 integrator */
    //SX sym_integrator = kmath::rk4_symbolic(X, U, dyn_func, dT);
    //Function RK4_INT = Function("RK4", {X,U,dT},{sym_integrator});

    /** CVODES integrator */
    /** @todo: make smarter initialisation of integrator */
    //double h = AlgoProps.sampling_time;
    //SX int_parameters = SX::vertcat({control, params});
    //SXDict ode = {{"x", state}, {"p", int_parameters}, {"ode", dynamics}};
    //Dict opts = {{"tf", h}};
    //Function CVODES_INT = integrator("CVODES_INT", "cvodes", ode, opts);

    /** assign class atributes */
    this->State = state;
    this->Control = control_cmd;
    this->Parameters = params;
    this->SymDynamics = dynamics;
    //this->SymIntegartor = sym_integrator;
    this->SymJacobian = d_jacobian;

    this->NumDynamics = dyn_func;
    this->NumSpecNongravForce = specNongravForce_func;
    this->NumJacobian = dyn_jac;

    /** return integrator function */
    /**
    if(AlgoProps.Integrator == IntType::CVODES)
        this->NumIntegrator = CVODES_INT;
    else
        this->NumIntegrator = RK4_INT;
        */
}


/** --------------------- */
/** Rigid Body Kinematics */
/** --------------------- */

RigidBodyKinematics::RigidBodyKinematics(const AlgorithmProperties &AlgoProps) {
    algo_props = AlgoProps;

    /** define system dynamics */
    SX r = SX::sym("r", 3);
    SX q = SX::sym("q", 4);
    SX vb = SX::sym("vb", 3);
    SX wb = SX::sym("wb", 3);
    SX u = SX::sym("u", 3);

    /** translation */
    SX qv = kmath::quat_multiply(q, SX::vertcat({0, vb}));
    SX qv_q = kmath::quat_multiply(qv, kmath::quat_inverse(q));
    SX rdot = qv_q(Slice(1, 4), 0);

    /** rotation / with norm correction*/
    double lambda = -10;
    SX qdot = 0.5 * kmath::quat_multiply(q, SX::vertcat({0, wb})) + 0.5 * lambda * q * (SX::dot(q, q) - 1);

    /** velocities */
    SX vb_dot = SX::zeros(3);
    SX wb_dot = SX::zeros(3);

    state = SX::vertcat({vb, wb, r, q});
    SX Dynamics = SX::vertcat({vb_dot, wb_dot, rdot, qdot});
    NumDynamics = Function("RB_Dynamics", {state}, {Dynamics});

    /** Jacobian */
    SX Jacobian = SX::jacobian(Dynamics, state);
    NumJacobian = Function("RB_Jacobian", {state, u}, {Jacobian});

    /** Integrators */
    /** CVODES */
    double h = algo_props.sampling_time;
    SXDict ode = {{"x",   state},
                  {"ode", Dynamics}};
    Dict opts = {{"tf", h}};
    Function CVODES_INT = integrator("CVODES_INT", "cvodes", ode, opts);
    NumIntegartor = CVODES_INT;
}

