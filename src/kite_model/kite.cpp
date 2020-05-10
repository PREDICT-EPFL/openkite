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


/* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder*/
template<typename W, typename GEN, typename DLO, typename DLA, typename AIL, typename ELV, typename RUD>
void KiteDynamics::getModel(GEN &g, GEN &rho,
                            W &windFrom, W &windSpeed,
                            GEN &b, GEN &c, GEN &AR, GEN &S,
                            GEN &Mass, GEN &Ixx, GEN &Iyy, GEN &Izz, GEN &Ixz,

                            GEN &e_oswald,
                            DLO &CD0,


                            DLO &CL0,
                            DLO &CLa,

                            DLO &Cm0,
                            DLO &Cma,


                            DLA &CYb,

                            GEN &Cl0,
                            DLA &Clb,

                            DLA &Cn0,
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

                            SX &v, SX &w, SX &r, SX &q,
                            SX &T, SX &dE, SX &dR, SX &dA,
                            SX &v_dot, SX &w_dot, SX &r_dot, SX &q_dot,
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

    /** Control variables **/
    T = SX::sym("T");   // Static propeller thrust along body frame x axis [N]
    dE = SX::sym("dE"); // Elevator deflection [positive causing negative pitch movement (nose down)] [rad]
    dR = SX::sym("dR"); // Rudder deflection [positive causing negative yaw movement (nose left)] [rad]
    dA = SX::sym("dA"); // Aileron deflection [positive causing negative roll movement (right wing up)] [rad]

    /** Aerodynamic (Wind) frame **/
    /* Wind velocity in body frame */
    SX g_vW = windSpeed * SX::vertcat({-cos(windFrom), -sin(windFrom), 0});
    SX b_vW = kmath::quat_transform(q_bg, g_vW);

    /* Apparent velocity in body frame, airspeed */
    SX b_va = v - b_vW;
    Va = SX::norm_2(b_va);

    /* Measured airspeed component (pitot tube orientation dependent) */
    SX r_sens = SX::vertcat({0.11, 0.22, -0.05});
    SX b_va_meas = b_va + SX::cross(r_sens, w);
    /* At fast body yawrates, the body rotation (thus pitot tube is faster than the CoG) is clearly negligible
     * under the effect of the pitot measurement direction being rotated out of the airflow (sideslip) */

    /* Pitot tube is oriented about 5 degress above body x axes */
    SX q_sens_b = kmath::T2quat(5.0 * M_PI / 180.0);
    SX sens_va = kmath::quat_transform(q_sens_b, b_va_meas);
    Va_pitot = sens_va(0);

    /* Aerodynamic angles (angle of attack, side slip angle */
    alpha = atan(b_va(2) / (b_va(0))); // + 1e-4));
    beta = asin(b_va(1) / (Va)); // + 1e-4));

    /** ---------------------------------------------------------- **/
    /** Aerodynamic Forces and Moments in aerodynamic (wind) frame **/
    /** ---------------------------------------------------------- **/
    SX dyn_press = 0.5 * rho * Va * Va;
    SX CL = CL0 + CLa * alpha + CLq * c / (2.0 * Va) * w(1) + CLde * dE;

    /** Forces in x, y, z directions: -Drag, Side force, -Lift **/
    SX LIFT = dyn_press * S * CL;

    SX DRAG = dyn_press * S * (CD0 + CL * CL / (pi * e_oswald * AR));

    SX SF = dyn_press * S * (CYb * beta +
                             b / (2.0 * Va) * (CYp * w(0) + CYr * w(2)) +
                             CYdr * dR);

    SX s_Faero = SX::vertcat({-DRAG, SF, -LIFT});

    /** Moments about x, y, z axes: L, M, N **/
    SX L = dyn_press * S * b * (Cl0 + Clb * beta +
                                b / (2.0 * Va) * (Clp * w(0) + Clr * w(2)) +
                                Clda * dA + Cldr * dR);

    SX M = dyn_press * S * c * (Cm0 + Cma * alpha +
                                c / (2.0 * Va) * Cmq +
                                Cmde * dE);

    SX N = dyn_press * S * b * (Cn0 + Cnb * beta +
                                b / (2.0 * Va) * (Cnp * w(0) + Cnr * w(2)) +
                                Cnda * dA + Cndr * dR);

    SX s_Maero = SX::vertcat({L, M, N});

    /** Aerodynamic Forces and Moments in body frame **/
    /* XFLR5 gives coefficients in stability frame */
    // To get from stability axis to body, rotate by aoa
    SX q_bs = kmath::T2quat(alpha);

    b_F_aero = kmath::quat_transform(q_bs, s_Faero);
    auto b_Maero = kmath::quat_transform(q_bs, s_Maero);

    /** ---------------------------------------- **/
    /** Gravitation, Thrust, Tether (body frame) **/
    /** ---------------------------------------- **/
    /** Gravitational acceleration **/
    SX b_g = kmath::quat_transform(q_bg, SX::vertcat({0, 0, g}));

    /** Propeller thrust **/
    const double p1 = -0.00083119;
    const double p2 = -0.014700129;
    b_F_thrust =
            T * (p1 * Va * Va + p2 * Va + 1.0) * SX::vertcat({1, 0, 0});  // T is the static thrust (at zero airspeed)

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
        SX decl = SX::atan(SX::norm_2(r(Slice(0, 2))) / -r(2));
        SX g_Wteth = 0.5 * (1 + cos(decl)) * tethLen * tethMassPerMeter * SX::vertcat({0, 0, g});
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
    q_dot = 0.5 * kmath::quat_multiply(q, SX::vertcat({0, w})) +    // q = q_gb, w = omega (body frame)
            0.5 * lambda * q * (SX::dot(q, q) - 1);                 // Quaternion norm stabilization term

    /** End of dynamics model ====================================================================================== **/
}

KiteDynamics::KiteDynamics(const KiteProperties &kiteProps, const AlgorithmProperties &AlgoProps,
                           const bool teth_ON, const bool controlsIncludeWind) {

    /** enviromental constants */
    double g = 9.80665; /** gravitational acceleration [m/s2] [WGS84] */
    double rho = 1.2985; /** standard atmospheric density [kg/m3] [standard Atmosphere 1976] */

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
    /**        Tether parameters       **/
    /** ------------------------------ **/
    //    double Lt = KiteProps.Tether.length;
    //    double Ks = KiteProps.Tether.Ks;
    //    double Kd = KiteProps.Tether.Kd;
    //    double rx = KiteProps.Tether.rx;
    //    double ry = KiteProps.Tether.ry;
    //    double rz = KiteProps.Tether.rz;


    SX v, w, r, q;
    SX T, dE, dR, dA;
    SX v_dot, w_dot, r_dot, q_dot;
    SX Va_pitot, Va, alpha, beta;
    SX b_F_aero, b_F_thrust, b_F_tether;

    SX control;

    if (controlsIncludeWind) {
        SX windFrom = SX::sym("windFrom");
        SX windSpeed = SX::sym("windSpeed");

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W,  GEN,   DLO,     DLA,   AIL,    ELV,    RUD*/
        getModel<SX, double, double, double, double, double, double>(
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

                v, w, r, q,
                T, dE, dR, dA,
                v_dot, w_dot, r_dot, q_dot,
                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, teth_ON, b_F_tether);

        control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

    } else {
        double windFrom = kiteProps.Wind.WindFrom;
        double windSpeed = kiteProps.Wind.WindSpeed;

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

                v, w, r, q,
                T, dE, dR, dA,
                v_dot, w_dot, r_dot, q_dot,
                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, teth_ON, b_F_tether);

        control = SX::vertcat({T, dE, dR, dA});
    }

    /** Complete dynamics of the Kite */
    auto state = SX::vertcat({v, w, r, q});
    auto dynamics = SX::vertcat({v_dot, w_dot, r_dot, q_dot});

    SX aero_values = SX::vertcat({Va, alpha, beta});

    Function dyn_func = Function("dynamics", {state, control}, {dynamics});
    Function airspeedMeas_func = Function("airspeed", {state}, {Va_pitot});
    Function aeroValues_func = Function("aero_out", {state}, {aero_values});
    Function specNongravForce_func = Function("spec_nongrav_force", {state, control},
                                              {(b_F_aero + b_F_thrust + b_F_tether) / Mass});
    Function specTethForce_func = Function("specTethForce", {state, control}, {b_F_tether / Mass});

    /** compute dynamics state Jacobian */
    SX d_jacobian = SX::jacobian(dynamics, state);
    Function dyn_jac = Function("dyn_jacobian", {state, control}, {d_jacobian});

    AeroDynamics = Function("Aero", {state, control}, {b_F_aero});
    /** define RK4 integrator scheme */
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control.size1());
    SX dT = SX::sym("dT");

    /** get symbolic expression for RK4 integrator */
    SX sym_integrator = kmath::rk4_symbolic(X, U, dyn_func, dT);
    Function RK4_INT = Function("RK4", {X, U, dT}, {sym_integrator});

    /** CVODES integrator */
    /** @todo: make smarter initialisation of integrator */
    double h = AlgoProps.sampling_time;
    SXDict ode = {{"x",   state},
                  {"p",   control},
                  {"ode", dynamics}};
    Dict opts = {{"tf", h}};
    Function CVODES_INT = integrator("CVODES_INT", "cvodes", ode, opts);

    /** assign class atributes */
    this->State = state;
    this->Control = control;
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
                           const kite_utils::IdentMode &identMode, const bool controlsIncludeWind) {

    /** enviromental constants */
    double g = 9.80665; /** gravitational acceleration [m/s2] [WGS84] */
    double rho = 1.2985; /** standard atmospheric density [kg/m3] [Standard Atmosphere 1976] */

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

    SX v, w, r, q;
    SX T, dE, dR, dA;
    SX v_dot, w_dot, r_dot, q_dot;
    SX Va_pitot, Va, alpha, beta;
    SX b_F_aero, b_F_thrust, b_F_tether;

    SX params;
    SX control;

    if (identMode == kite_utils::IdentMode::LONGITUDINAL) {
        /** LONGITUDINAL IDENTIFICATION PARAMETERS ------------------------------------------------------------------ */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
        double e_oswald = kiteProps.Aerodynamics.e_oswald;
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

        params = SX::vertcat({CD0,

                              CL0,
                              CLa,

                              Cm0,
                              Cma,

                              CLq,
                              Cmq,

                              CLde,
                              Cmde
                             }); // 9 longitudinal parameters

        if (controlsIncludeWind) {
            SX windFrom = SX::sym("windFrom");
            SX windSpeed = SX::sym("windSpeed");

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W, GEN,   DLO, DLA,   AIL,    ELV, RUD*/
            getModel<SX, double, SX, double, double, SX, double>(
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

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

        } else {
            double windFrom = kiteProps.Wind.WindFrom;
            double windSpeed = kiteProps.Wind.WindSpeed;

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

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA});
        }

    } else if (identMode == kite_utils::IdentMode::LATERAL) {
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

        SX Cn0 = SX::sym("Cn0");
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

        params = SX::vertcat({params,
                              CYb,

                              Clb,

                              Cn0,
                              Cnb,

                              CYp,
                              Clp,
                              Cnp,

                              CYr,
                              Clr,
                              Cnr,

                              Clda,
                              Cnda
                             }); // 12 lateral parameters (Aileron control only)

        if (controlsIncludeWind) {
            SX windFrom = SX::sym("windFrom");
            SX windSpeed = SX::sym("windSpeed");

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W, GEN,    DLO,   DLA, AIL, ELV,   RUD*/
            getModel<SX, double, double, SX, SX, double, double>(
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

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

        } else {
            double windFrom = kiteProps.Wind.WindFrom;
            double windSpeed = kiteProps.Wind.WindSpeed;

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

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA});
        }
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

        SX Cn0 = SX::sym("Cn0");
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

        params = SX::vertcat({params,
                              CYb,

                              Clb,

                              Cn0,
                              Cnb,

                              CYp,
                              Clp,
                              Cnp,

                              CYr,
                              Clr,
                              Cnr,

                              CYdr,
                              Cldr,
                              Cndr
                             }); // 13 lateral parameters (Rudder control only)

        if (controlsIncludeWind) {
            SX windFrom = SX::sym("windFrom");
            SX windSpeed = SX::sym("windSpeed");

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W, GEN,    DLO,   DLA, AIL,    ELV,   RUD*/
            getModel<SX, double, double, SX, double, double, SX>(
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

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

        } else {
            double windFrom = kiteProps.Wind.WindFrom;
            double windSpeed = kiteProps.Wind.WindSpeed;

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

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA});
        }


    } else if (identMode == kite_utils::IdentMode::COMPLETE) {
        /** COMPLETE IDENTIFICATION PARAMETERS ------------------------------------------------------------------ */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
        double e_oswald = kiteProps.Aerodynamics.e_oswald;
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

        SX Cn0 = SX::sym("Cn0");
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

        params = SX::vertcat({CD0,

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
                              Cndr
                             }); // 25 parameters


        if (controlsIncludeWind) {
            SX windFrom = SX::sym("windFrom");
            SX windSpeed = SX::sym("windSpeed");

            /* Info: <General, Lon, Lat> Types for parameter groups, defined in function call */
            getModel<SX, double, SX, SX>(g, rho,
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

                                         v, w, r, q,
                                         T, dE, dR, dA,
                                         v_dot, w_dot, r_dot, q_dot,
                                         Va_pitot, Va, alpha, beta,
                                         b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

        } else {
            double windFrom = kiteProps.Wind.WindFrom;
            double windSpeed = kiteProps.Wind.WindSpeed;

            /* Info: <General, Lon, Lat> Types for parameter groups, defined in function call */
            getModel<double, double, SX, SX>(g, rho,
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

                                             v, w, r, q,
                                             T, dE, dR, dA,
                                             v_dot, w_dot, r_dot, q_dot,
                                             Va_pitot, Va, alpha, beta,
                                             b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA});
        }
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
    auto state = SX::vertcat({v, w, r, q});
    auto dynamics = SX::vertcat({v_dot, w_dot, r_dot, q_dot});

    Function dyn_func = Function("dynamics", {state, control, params}, {dynamics});
    Function specNongravForce_func = Function("spec_nongrav_force", {state, control}, {(b_F_aero + b_F_thrust) / Mass});

    /** compute dynamics state Jacobian */
    SX d_jacobian = SX::jacobian(dynamics, state);
    Function dyn_jac = Function("dyn_jacobian", {state, control, params}, {d_jacobian});

    /** define RK4 integrator scheme */
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control.size1());
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
    this->Control = control;
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


/* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder*/
template<typename W, typename GEN, typename DLO, typename DLA, typename AIL, typename ELV, typename RUD>
void MinimalKiteDynamics::getMinimalModel(GEN &g, GEN &rho,
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

        //GEN &Cl0,
                                          DLA &Clb,

        //GEN &Cn0,
                                          DLA &Cnb,


        //DLO &CLq,
        //DLO &Cmq,

        //GEN &CYp,
                                          DLA &Clp,
        //DLA &Cnp,

        //DLA &CYr,
        //DLA &Clr,
                                          DLA &Cnr,


        //ELV &CLde,
                                          ELV &Cmde,

                                          AIL &Clda,
        //GEN &Cnda,

        //RUD &CYdr,
        //RUD &Cldr,
                                          RUD &Cndr,

                                          SX &v, SX &w, SX &r, SX &q,
                                          SX &T, SX &dE, SX &dR, SX &dA,
                                          SX &v_dot, SX &w_dot, SX &r_dot, SX &q_dot,
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

    /** Control variables **/
    T = SX::sym("T");   // Static propeller thrust along body frame x axis [N]
    dE = SX::sym("dE"); // Elevator deflection [positive causing negative pitch movement (nose down)] [rad]
    dR = SX::sym("dR"); // Rudder deflection [positive causing negative yaw movement (nose left)] [rad]
    dA = SX::sym("dA"); // Aileron deflection [positive causing negative roll movement (right wing up)] [rad]

    /** Aerodynamic (Wind) frame **/
    /* Wind velocity in body frame */
    SX g_vW = windSpeed * SX::vertcat({-cos(windFrom), -sin(windFrom), 0});
    SX b_vW = kmath::quat_transform(q_bg, g_vW);

    /* Apparent velocity in body frame, airspeed */
    SX b_va = v - b_vW;
    Va = SX::norm_2(b_va);

    /* Measured airspeed component (pitot tube orientation dependent) */
    /* //    SX r_sens = SX::vertcat({0.11, 0.21, -0.05});
     * //    SX b_va_meas = b_va + SX::cross(r_sens, w);
     * // At fast body yawrates. the body rotation (thus pitot tube is faster than the CoG) is clearly negligible
     * // under the effect of the pitot measurement direction being rotated out of the airflow (sideslip)
     * */

    /* //    SX q_sens_b = kmath::T2quat(5.0 * M_PI / 180.0);
     * //    SX sens_va = kmath::quat_transform(q_sens_b, b_va_meas);
     * //    Va_pitot = sens_va(0);
     * // Also, the 5-degree angle of the pitot makes 0.5 percent error. So we take the body x-direction as
     * airspeed measurement, which decreases at any sidelip angle.
     * */
    Va_pitot = b_va(0);

    /* Aerodynamic angles (angle of attack, side slip angle */
    SX aoa = atan(b_va(2) / (b_va(0) + 1e-4));
    SX ss = asin(b_va(1) / (Va + 1e-4));

    /** ---------------------------------------------------------- **/
    /** Aerodynamic Forces and Moments in stability frame **/
    /** ---------------------------------------------------------- **/
    SX dyn_press = 0.5 * rho * Va * Va;
    SX CL = CL0 + CLa * aoa;

    /** Forces in x, y, z directions: -Drag, Side force, -Lift **/
    SX LIFT = dyn_press * S * CL;

    SX DRAG = dyn_press * S * (CD0 + CL * CL / (pi * e_oswald * AR));

    SX SF = dyn_press * S * CYb * ss;

    SX s_Faero = SX::vertcat({-DRAG, SF, -LIFT});

    /** Moments about x, y, z axes: L, M, N **/
    SX L = dyn_press * S * b * (Clb * ss +
                                b / (2.0 * Va) * Clp * w(0) +
                                Clda * dA);

    SX M = dyn_press * S * c * (Cm0 + Cma * aoa +
                                Cmde * dE);

    SX N = dyn_press * S * b * (Cnb * ss +
                                b / (2.0 * Va) * Cnr * w(2) +
                                Cndr * dR);

    SX s_Maero = SX::vertcat({L, M, N});

    /** Aerodynamic Forces and Moments in body frame **/
    /* XFLR5 gives coefficients in stability frame */
    // To get from stability axis to body, rotate by aoa
    SX q_bs = kmath::T2quat(aoa);

    b_F_aero = kmath::quat_transform(q_bs, s_Faero);
    auto b_Maero = kmath::quat_transform(q_bs, s_Maero);

    /** ---------------------------------------- **/
    /** Gravitation, Thrust, Tether (body frame) **/
    /** ---------------------------------------- **/
    /** Gravitational acceleration **/
    SX b_g = kmath::quat_transform(q_bg, SX::vertcat({0, 0, g}));

    /** Propeller thrust **/
    const double p1 = -0.00083119;
    const double p2 = -0.014700129;
    b_F_thrust =
            T * (p1 * Va * Va + p2 * Va + 1.0) * SX::vertcat({1, 0, 0});  // T is the static thrust (at zero airspeed)

    /** Tether force and moment **/
    b_F_tether = SX::vertcat({0, 0, 0});

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
    q_dot = 0.5 * kmath::quat_multiply(q, SX::vertcat({0, w})) +    // q = q_gb, w = omega (body frame)
            0.5 * lambda * q * (SX::dot(q, q) - 1);                 // Quaternion norm stabilization term

    /** End of dynamics model ====================================================================================== **/
}

MinimalKiteDynamics::MinimalKiteDynamics(const KiteProperties &kiteProps, const AlgorithmProperties &AlgoProps,
                                         const bool teth_ON, const bool controlsIncludeWind) {

    /** enviromental constants */
    double g = 9.80665; /** gravitational acceleration [m/s2] [WGS84] */
    double rho = 1.2985; /** standard atmospheric density [kg/m3] [standard Atmosphere 1976] */

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

    double Clb = kiteProps.Aerodynamics.Clb;

    double Cnb = kiteProps.Aerodynamics.Cnb;

    /* Pitchrate */
//    double Cmq = kiteProps.Aerodynamics.Cmq;

    /* Rollrate */
    double Clp = kiteProps.Aerodynamics.Clp;
//    double Cnp = kiteProps.Aerodynamics.Cnp;

    /* Yawrate */
//    double Clr = kiteProps.Aerodynamics.Clr;
    double Cnr = kiteProps.Aerodynamics.Cnr;


    /** ------------------------------ **/
    /** Aerodynamic effects of control **/
    /** ------------------------------ **/
    /* Elevator */
    double Cmde = kiteProps.Aerodynamics.Cmde;

    /* Ailerons */
    double Clda = kiteProps.Aerodynamics.Clda;

    /* Rudder */
    double Cndr = kiteProps.Aerodynamics.Cndr;

    /** ------------------------------ **/
    /**        Tether parameters       **/
    /** ------------------------------ **/
    //    double Lt = KiteProps.Tether.length;
    //    double Ks = KiteProps.Tether.Ks;
    //    double Kd = KiteProps.Tether.Kd;
    //    double rx = KiteProps.Tether.rx;
    //    double ry = KiteProps.Tether.ry;
    //    double rz = KiteProps.Tether.rz;


    SX v, w, r, q;
    SX T, dE, dR, dA;
    SX v_dot, w_dot, r_dot, q_dot;
    SX Va_pitot, Va, alpha, beta;
    SX b_F_aero, b_F_thrust, b_F_tether;

    SX control;

    if (controlsIncludeWind) {
        SX windFrom = SX::sym("windFrom");
        SX windSpeed = SX::sym("windSpeed");

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W,  GEN,   DLO,     DLA,   AIL,    ELV,    RUD*/
        getMinimalModel<SX, double, double, double, double, double, double>(
                g, rho,
                windFrom, windSpeed,
                b, c, AR, S,
                Mass, Ixx, Iyy, Izz, Ixz,

                e_oswald, CD0,

                CL0, CLa,
                Cm0, Cma,

                CYb,
                Clb,
                Cnb,

                Clp,

                Cnr,

                Cmde,
                Clda,
                Cndr,

                v, w, r, q,
                T, dE, dR, dA,
                v_dot, w_dot, r_dot, q_dot,
                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, teth_ON, b_F_tether);

        control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

    } else {
        double windFrom = kiteProps.Wind.WindFrom;
        double windSpeed = kiteProps.Wind.WindSpeed;

        /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
         *        W,     GEN,   DLO,     DLA,   AIL,    ELV,    RUD*/
        getMinimalModel<double, double, double, double, double, double, double>(
                g, rho,
                windFrom, windSpeed,
                b, c, AR, S,
                Mass, Ixx, Iyy, Izz, Ixz,

                e_oswald, CD0,

                CL0, CLa,
                Cm0, Cma,

                CYb,
                Clb,
                Cnb,

                Clp,

                Cnr,

                Cmde,
                Clda,
                Cndr,

                v, w, r, q,
                T, dE, dR, dA,
                v_dot, w_dot, r_dot, q_dot,
                Va_pitot, Va, alpha, beta,
                b_F_aero, b_F_thrust, teth_ON, b_F_tether);

        control = SX::vertcat({T, dE, dR, dA});
    }

    /** Complete dynamics of the Kite */
    auto state = SX::vertcat({v, w, r, q});
    auto dynamics = SX::vertcat({v_dot, w_dot, r_dot, q_dot});

    SX aero_values = SX::vertcat({Va, alpha, beta});

    Function dyn_func = Function("dynamics", {state, control}, {dynamics});
    Function airspeed_meas_func = Function("airspeed", {state, control}, {Va_pitot});
    Function aero_values_func = Function("aero_out", {state}, {aero_values});
    Function specNongravForce_func = Function("spec_nongrav_force", {state, control},
                                              {(b_F_aero + b_F_thrust + b_F_tether) / Mass});
    Function specTethForce_func = Function("specTethForce", {state, control}, {b_F_tether / Mass});

    /** compute dynamics state Jacobian */
    SX d_jacobian = SX::jacobian(dynamics, state);
    Function dyn_jac = Function("dyn_jacobian", {state, control}, {d_jacobian});

    AeroDynamics = Function("Aero", {state, control}, {b_F_aero});
    /** define RK4 integrator scheme */
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control.size1());
    SX dT = SX::sym("dT");

    /** get symbolic expression for RK4 integrator */
    SX sym_integrator = kmath::rk4_symbolic(X, U, dyn_func, dT);
    Function RK4_INT = Function("RK4", {X, U, dT}, {sym_integrator});

    /** CVODES integrator */
    /** @todo: make smarter initialisation of integrator */
    double h = AlgoProps.sampling_time;
    SXDict ode = {{"x",   state},
                  {"p",   control},
                  {"ode", dynamics}};
    Dict opts = {{"tf", h}};
    Function CVODES_INT = integrator("CVODES_INT", "cvodes", ode, opts);

    /** assign class atributes */
    this->State = state;
    this->Control = control;
    this->SymDynamics = dynamics;
    this->SymIntegartor = sym_integrator;
    this->SymJacobian = d_jacobian;

    this->NumDynamics = dyn_func;
    this->NumAirspeedMeasured = airspeed_meas_func;
    this->NumAeroValues = aero_values_func;
    this->NumSpecNongravForce = specNongravForce_func;
    this->NumSpecTethForce = specTethForce_func;
    this->NumJacobian = dyn_jac;

    /** return integrator function */
    if (AlgoProps.Integrator == IntType::CVODES)
        this->NumIntegrator = CVODES_INT;
    else
        this->NumIntegrator = RK4_INT;

}

MinimalKiteDynamics::MinimalKiteDynamics(const KiteProperties &kiteProps, const AlgorithmProperties &AlgoProps,
                                         const kite_utils::IdentMode &identMode, const bool controlsIncludeWind) {

    /** enviromental constants */
    double g = 9.80665; /** gravitational acceleration [m/s2] [WGS84] */
    double rho = 1.2985; /** standard atmospheric density [kg/m3] [Standard Atmosphere 1976] */

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

    SX v, w, r, q;
    SX T, dE, dR, dA;
    SX v_dot, w_dot, r_dot, q_dot;
    SX Va_pitot, Va, alpha, beta;
    SX b_F_aero, b_F_thrust, b_F_tether;

    SX params;
    SX control;

    if (identMode == kite_utils::IdentMode::LONGITUDINAL) {
        /** LONGITUDINAL IDENTIFICATION PARAMETERS ------------------------------------------------------------------ */

        /** ------------------------------- **/
        /** Aerodynamic parameters          **/
        /** ------------------------------- **/
        SX e_oswald = SX::sym("e_oswald");
        SX CD0 = SX::sym("CD0");

        /* AOA */
        SX CL0 = SX::sym("CL0");
        SX CLa = SX::sym("CLa");

        SX Cm0 = SX::sym("Cm0");
        SX Cma = SX::sym("Cma");

        /* Sideslip */
        double CYb = kiteProps.Aerodynamics.CYb;

//        double Cl0 = kiteProps.Aerodynamics.Cl0;
        double Clb = kiteProps.Aerodynamics.Clb;

//        double Cn0 = kiteProps.Aerodynamics.Cn0;
        double Cnb = kiteProps.Aerodynamics.Cnb;


        /* Pitchrate */
//        SX CLq = SX::sym("CLq");
//        SX Cmq = SX::sym("Cmq");

        /* Rollrate */
//        double CYp = kiteProps.Aerodynamics.CYp;
        double Clp = kiteProps.Aerodynamics.Clp;
//        double Cnp = kiteProps.Aerodynamics.Cnp;

        /* Yawrate */
//        double CYr = kiteProps.Aerodynamics.CYr;
//        double Clr = kiteProps.Aerodynamics.Clr;
        double Cnr = kiteProps.Aerodynamics.Cnr;


        /** ------------------------------ **/
        /** Aerodynamic effects of control **/
        /** ------------------------------ **/
        /* Elevator */
//        SX CLde = SX::sym("CLde");
        SX Cmde = SX::sym("Cmde");

        /* Ailerons */
        double Clda = kiteProps.Aerodynamics.Clda;
//        double Cnda = kiteProps.Aerodynamics.Cnda;

        /* Rudder */
//        double CYdr = kiteProps.Aerodynamics.CYdr;
//        double Cldr = kiteProps.Aerodynamics.Cldr;
        double Cndr = kiteProps.Aerodynamics.Cndr;

        params = SX::vertcat({e_oswald,
                              CD0,

                              CL0,
                              CLa,

                              Cm0,
                              Cma,

//                              CLq,
//                              Cmq,
//
//                              CLde,
                              Cmde
                             }); // 7 longitudinal parameters

        if (controlsIncludeWind) {
            SX windFrom = SX::sym("windFrom");
            SX windSpeed = SX::sym("windSpeed");

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W, GEN,   DLO, DLA,   AIL,    ELV, RUD*/
            getMinimalModel<SX, double, SX, double, double, SX, double>(
                    g, rho,
                    windFrom, windSpeed,
                    b, c, AR, S,
                    Mass, Ixx, Iyy, Izz, Ixz,

                    e_oswald, CD0,

                    CL0, CLa,
                    Cm0, Cma,

                    CYb,
                    Clb,
                    Cnb,

                    Clp,

                    Cnr,

                    Cmde,
                    Clda,
                    Cndr,

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

        } else {
            double windFrom = kiteProps.Wind.WindFrom;
            double windSpeed = kiteProps.Wind.WindSpeed;

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W, GEN,   DLO, DLA,   AIL,    ELV, RUD*/
            getMinimalModel<double, double, SX, double, double, SX, double>(
                    g, rho,
                    windFrom, windSpeed,
                    b, c, AR, S,
                    Mass, Ixx, Iyy, Izz, Ixz,

                    e_oswald, CD0,

                    CL0, CLa,
                    Cm0, Cma,

                    CYb,
                    Clb,
                    Cnb,

                    Clp,

                    Cnr,

                    Cmde,
                    Clda,
                    Cndr,

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA});
        }

    } else if (identMode == kite_utils::IdentMode::LATERAL) {
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

//        double Cl0 = kiteProps.Aerodynamics.Cl0;
        SX Clb = SX::sym("Clb");

//        double Cn0 = kiteProps.Aerodynamics.Cn0;
        SX Cnb = SX::sym("Cnb");


        /* Pitchrate */
//        double CLq = kiteProps.Aerodynamics.CLq;
//        double Cmq = kiteProps.Aerodynamics.Cmq;

        /* Rollrate */
        //SX CYp = SX::sym("CYp");
//        double CYp = kiteProps.Aerodynamics.CYp;
        SX Clp = SX::sym("Clp");
//        SX Cnp = SX::sym("Cnp");

        /* Yawrate */
//        SX CYr = SX::sym("CYr");
//        SX Clr = SX::sym("Clr");
        SX Cnr = SX::sym("Cnr");


        /** ------------------------------ **/
        /** Aerodynamic effects of control **/
        /** ------------------------------ **/
        /* Elevator */
//        double CLde = kiteProps.Aerodynamics.CLde;
        double Cmde = kiteProps.Aerodynamics.Cmde;

        /* Ailerons */
        SX Clda = SX::sym("Clda");
        //SX Cnda = SX::sym("Cnda");
//        double Cnda = kiteProps.Aerodynamics.Cnda;

        /* Rudder */
//        double CYdr = kiteProps.Aerodynamics.CYdr;
//        double Cldr = kiteProps.Aerodynamics.Cldr;
        double Cndr = kiteProps.Aerodynamics.Cndr;

        params = SX::vertcat({params,
                              CYb,

                              Clb,

                              Cnb,

                                     //CYp,
                              Clp,
                                     //Cnp,

                                     //CYr,
                                     //Clr,
                              Cnr,

                              Clda,
                                     //Cnda
                             }); // 6 lateral parameters (Aileron control only)

        if (controlsIncludeWind) {
            SX windFrom = SX::sym("windFrom");
            SX windSpeed = SX::sym("windSpeed");

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W, GEN,    DLO,   DLA, AIL, ELV,   RUD*/
            getMinimalModel<SX, double, double, SX, SX, double, double>(
                    g, rho,
                    windFrom, windSpeed,
                    b, c, AR, S,
                    Mass, Ixx, Iyy, Izz, Ixz,

                    e_oswald, CD0,

                    CL0, CLa,
                    Cm0, Cma,

                    CYb,
                    Clb,
                    Cnb,

                    Clp,

                    Cnr,

                    Cmde,
                    Clda,
                    Cndr,

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

        } else {
            double windFrom = kiteProps.Wind.WindFrom;
            double windSpeed = kiteProps.Wind.WindSpeed;

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W,     GEN,    DLO,   DLA, AIL, ELV,   RUD*/
            getMinimalModel<double, double, double, SX, SX, double, double>(
                    g, rho,
                    windFrom, windSpeed,
                    b, c, AR, S,
                    Mass, Ixx, Iyy, Izz, Ixz,

                    e_oswald, CD0,

                    CL0, CLa,
                    Cm0, Cma,

                    CYb,
                    Clb,
                    Cnb,

                    Clp,

                    Cnr,

                    Cmde,
                    Clda,
                    Cndr,

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA});
        }
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

//        double Cl0 = kiteProps.Aerodynamics.Cl0;
        SX Clb = SX::sym("Clb");

//        double Cn0 = kiteProps.Aerodynamics.Cn0;
        SX Cnb = SX::sym("Cnb");


        /* Pitchrate */
//        double CLq = kiteProps.Aerodynamics.CLq;
//        double Cmq = kiteProps.Aerodynamics.Cmq;

        /* Rollrate */
        //SX CYp = SX::sym("CYp");
//        double CYp = kiteProps.Aerodynamics.CYp;
        SX Clp = SX::sym("Clp");
//        SX Cnp = SX::sym("Cnp");

        /* Yawrate */
//        SX CYr = SX::sym("CYr");
//        SX Clr = SX::sym("Clr");
        SX Cnr = SX::sym("Cnr");


        /** ------------------------------ **/
        /** Aerodynamic effects of control **/
        /** ------------------------------ **/
        /* Elevator */
//        double CLde = kiteProps.Aerodynamics.CLde;
        double Cmde = kiteProps.Aerodynamics.Cmde;

        /* Ailerons */
        double Clda = kiteProps.Aerodynamics.Clda;
//        double Cnda = kiteProps.Aerodynamics.Cnda;

        /* Rudder */
//        SX CYdr = SX::sym("CYdr");
//        SX Cldr = SX::sym("Cldr");
        SX Cndr = SX::sym("Cndr");

        params = SX::vertcat({params,
                              CYb,

                              Clb,

                              Cnb,

                                     //CYp,
                              Clp,
                                     //Cnp,

                                     //CYr,
                                     //Clr,
                              Cnr,

                                     //CYdr,
                                     //Cldr,
                              Cndr
                             }); // 6 lateral parameters (Rudder control only)

        if (controlsIncludeWind) {
            SX windFrom = SX::sym("windFrom");
            SX windSpeed = SX::sym("windSpeed");

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W, GEN,    DLO,   DLA, AIL,    ELV,   RUD*/
            getMinimalModel<SX, double, double, SX, double, double, SX>(
                    g, rho,
                    windFrom, windSpeed,
                    b, c, AR, S,
                    Mass, Ixx, Iyy, Izz, Ixz,

                    e_oswald, CD0,

                    CL0, CLa,
                    Cm0, Cma,

                    CYb,
                    Clb,
                    Cnb,

                    Clp,

                    Cnr,

                    Cmde,
                    Clda,
                    Cndr,

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});

        } else {
            double windFrom = kiteProps.Wind.WindFrom;
            double windSpeed = kiteProps.Wind.WindSpeed;

            /* Wind, GENeral, Dynamic LOngitudinal, Dynamic LAteral, AILeron, ELeVator, RUDder
             *        W,     GEN,    DLO,   DLA, AIL,    ELV,   RUD*/
            getMinimalModel<double, double, double, SX, double, double, SX>(
                    g, rho,
                    windFrom, windSpeed,
                    b, c, AR, S,
                    Mass, Ixx, Iyy, Izz, Ixz,

                    e_oswald, CD0,

                    CL0, CLa,
                    Cm0, Cma,

                    CYb,
                    Clb,
                    Cnb,

                    Clp,

                    Cnr,

                    Cmde,
                    Clda,
                    Cndr,

                    v, w, r, q,
                    T, dE, dR, dA,
                    v_dot, w_dot, r_dot, q_dot,
                    Va_pitot, Va, alpha, beta,
                    b_F_aero, b_F_thrust, false, b_F_tether);

            control = SX::vertcat({T, dE, dR, dA});
        }


    }

//    else if (identMode == COMPLETE) {
//        /** COMPLETE IDENTIFICATION PARAMETERS ------------------------------------------------------------------ */
//
//        /** ------------------------------- **/
//        /** Aerodynamic parameters          **/
//        /** ------------------------------- **/
//        double e_oswald = kiteProps.Aerodynamics.e_oswald;
//        SX CD0 = SX::sym("CD0");
//
//        /* AOA */
//        SX CL0 = SX::sym("CL0");
//        SX CLa = SX::sym("CLa");
//
//        SX Cm0 = SX::sym("Cm0");
//        SX Cma = SX::sym("Cma");
//
//        /* Sideslip */
//        SX CYb = SX::sym("CYb");
//
//        double Cl0 = kiteProps.Aerodynamics.Cl0;
//        SX Clb = SX::sym("Clb");
//
//        double Cn0 = kiteProps.Aerodynamics.Cn0;
//        SX Cnb = SX::sym("Cnb");
//
//
//        /* Pitchrate */
//        SX CLq = SX::sym("CLq");
//        SX Cmq = SX::sym("Cmq");
//
//        /* Rollrate */
//        SX CYp = SX::sym("CYp");
//        SX Clp = SX::sym("Clp");
//        SX Cnp = SX::sym("Cnp");
//
//        /* Yawrate */
//        SX CYr = SX::sym("CYr");
//        SX Clr = SX::sym("Clr");
//        SX Cnr = SX::sym("Cnr");
//
//
//        /** ------------------------------ **/
//        /** Aerodynamic effects of control **/
//        /** ------------------------------ **/
//        /* Elevator */
//        SX CLde = SX::sym("CLde");
//        SX Cmde = SX::sym("Cmde");
//
//        /* Ailerons */
//        SX Clda = SX::sym("Clda");
//        SX Cnda = SX::sym("Cnda");
//
//        /* Rudder */
//        double CYdr = kiteProps.Aerodynamics.CYdr;
//        double Cldr = kiteProps.Aerodynamics.Cldr;
//        double Cndr = kiteProps.Aerodynamics.Cndr;
//
//        params = SX::vertcat({CD0,
//
//                              CL0,
//                              CLa,
//
//                              Cm0,
//                              Cma,
//
//
//                              CYb,
//
//                              Clb,
//
//                              Cnb,
//
//
//                              CLq,
//                              Cmq,
//
//                              CYp,
//                              Clp,
//                              Cnp,
//
//                              CYr,
//                              Clr,
//                              Cnr,
//
//
//                              CLde,
//                              Cmde,
//
//                              Clda,
//                              Cnda
//                             }); // X parameters
//
//
//        if (controlsIncludeWind) {
//            SX windFrom = SX::sym("windFrom");
//            SX windSpeed = SX::sym("windSpeed");
//
//            /* Info: <General, Lon, Lat> Types for parameter groups, defined in function call */
//            getModel<SX, double, SX, SX>(g, rho,
//                                         windFrom, windSpeed,
//                                         b, c, AR, S,
//                                         Mass, Ixx, Iyy, Izz, Ixz,
//
//                                         e_oswald, CD0,
//
//                                         CL0, CLa,
//                                         Cm0, Cma,
//
//                                         CYb,
//                                         Cl0, Clb,
//                                         Cn0, Cnb,
//
//                                         CLq, Cmq,
//                                         CYp, Clp, Cnp,
//                                         CYr, Clr, Cnr,
//
//                                         CLde, Cmde,
//                                         Clda, Cnda,
//                                         CYdr, Cldr, Cndr,
//
//                                         v, w, r, q,
//                                         T, dE, dR, dA,
//                                         v_dot, w_dot, r_dot, q_dot,
//                                         Va_pitot, Va, alpha, beta,
//                                         b_F_aero, b_F_thrust, false, b_F_tether);
//
//            control = SX::vertcat({T, dE, dR, dA, windFrom, windSpeed});
//
//        } else {
//            double windFrom = kiteProps.Wind.WindFrom;
//            double windSpeed = kiteProps.Wind.WindSpeed;
//
//            /* Info: <General, Lon, Lat> Types for parameter groups, defined in function call */
//            getModel<double, double, SX, SX>(g, rho,
//                                             windFrom, windSpeed,
//                                             b, c, AR, S,
//                                             Mass, Ixx, Iyy, Izz, Ixz,
//
//                                             e_oswald, CD0,
//
//                                             CL0, CLa,
//                                             Cm0, Cma,
//
//                                             CYb,
//                                             Cl0, Clb,
//                                             Cn0, Cnb,
//
//                                             CLq, Cmq,
//                                             CYp, Clp, Cnp,
//                                             CYr, Clr, Cnr,
//
//                                             CLde, Cmde,
//                                             Clda, Cnda,
//                                             CYdr, Cldr, Cndr,
//
//                                             v, w, r, q,
//                                             T, dE, dR, dA,
//                                             v_dot, w_dot, r_dot, q_dot,
//                                             Va_pitot, Va, alpha, beta,
////                                         b_F_aero, b_F_thrust, false, b_F_tether);
//
//            control = SX::vertcat({T, dE, dR, dA});
//        }
//    }

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
    auto state = SX::vertcat({v, w, r, q});
    auto dynamics = SX::vertcat({v_dot, w_dot, r_dot, q_dot});

    Function dyn_func = Function("dynamics", {state, control, params}, {dynamics});
    Function specNongravForce_func = Function("spec_nongrav_force", {state, control}, {(b_F_aero + b_F_thrust) / Mass});

    /** compute dynamics state Jacobian */
    SX d_jacobian = SX::jacobian(dynamics, state);
    Function dyn_jac = Function("dyn_jacobian", {state, control, params}, {d_jacobian});

    /** define RK4 integrator scheme */
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control.size1());
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
    this->Control = control;
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

