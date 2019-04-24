#include "mobile_robot.h"

using namespace casadi;


casadi::SX reg_sign(const casadi::SX &x)
{
    return ((1 + 2 * sign(x)) / 2);
}


MobileRobot::MobileRobot(const MobileRobotProperties &props)
{
    /** model constants [DevBot2]*/
    double g    = 9.81;    // gravity acceleration

    double m    = props.mass;
    double b    = props.b;
    double a    = props.a;
    double Iz   = props.Iz;
    double r_f  = props.r_f;
    double r_r  = props.r_r;
    double fwIz = props.fwIz;
    double rwIz = props.rwIz;
    double Cy_r = props.Cy_r;
    double Cx_r = props.Cx_r;
    double Cx_f = props.Cx_f;
    double Cy_f = props.Cy_f;
    double mu   = props.mu;
    double rollResist = props.rollResist;
    double rollResistSpeed = props.rollResistSpeed;
    double Cxx = props.Cxx;
    TIRE_MODEL tire = LINEAR;//props.tire_model;



    /** state definition */
    SX vx    = SX::sym("vx");
    SX vy    = SX::sym("vy");
    SX omega = SX::sym("omega");

    SX x     = SX::sym("x");
    SX y     = SX::sym("y");
    SX theta = SX::sym("theta");

    SX owf   = SX::sym("owf");
    SX owr   = SX::sym("owr");
    state = SX::vertcat({vx, vy, omega, x, y, theta, owf, owr});

    /** controls definition */
    SX phi       = SX::sym("phi");
    SX FXr_req   = SX::sym("Fxr_req");
    SX FXf_req   = SX::sym("FXf_req");
    SX control   = SX::vertcat({phi, FXf_req, FXr_req});

    /** tire forces */
    SX FXf, FXr, FYf, FYr, Fdrag;

    /** vertical force acting on the wheels */
    SX FZr, FZf;

    /** REAR wheel force */
    SX v_wr = r_r * owr;
    SX vx_r = vx;
    SX vy_r = vy - omega * b;
    FZr = m * g * a / (a + b);

    if(tire == BRUSH)
    {
        SX tan_ar = tan( vy_r / (vx + reg_sign(vx) * 1e-3));   //regularization
        SX sigma_yr = -Cy_r * tan_ar * (vx / (norm_1(v_wr) + 1e-3));   // "regularization"
        SX sigma_xr = Cx_r * (v_wr - vx_r) / (norm_1(v_wr) + 1e-3);     //question Vsx sign
        SX gamma_r = norm_2(SX::vertcat({sigma_xr + reg_sign(sigma_xr) * 10, sigma_yr + reg_sign(sigma_yr) * 10}));

        /** Brush model of the tire */
        SX gamma_s = 3 * mu * FZr;
        SX f2 = mu * FZr;
        SX f1 = gamma_r - pow(gamma_r,2) / gamma_s + (1 / 3.0) * pow(gamma_r, 3) / pow(gamma_s, 2);

        SX Fs_r = (1 - heaviside(gamma_r - gamma_s + 1e-5)) * f1 + heaviside(gamma_r - gamma_s + 1e-5) * f2;
        FXr = (sigma_xr / (gamma_r)) * Fs_r;
        FYr = (sigma_yr / (gamma_r)) * Fs_r;
    }
    else if(tire == LINEAR)
    {
        /** longitudinal slip */
        SX k_r = (v_wr - vx_r) / (vx_r);     // should regularize here

        FXr = Cx_r * k_r;
        FYr = -Cy_r * atan2(vy_r, abs(vx_r));                    //atan(tan_ar);
    }

    /** FRONT wheel force */
    SX v_wf = r_f * owf;                                         // "regularization"
    SX vx_f = cos(-phi) * vx - sin(-phi) * (vy + omega * a);
    SX vy_f = sin(-phi) * vx + cos(-phi) * (vy + omega * a);
    FZf = m * g * b / (a + b);


    if(tire == BRUSH)
    {
        SX tan_af = tan( vy_f / (vx_f + reg_sign(vx_f) * 1e-3));
        SX sigma_yf = -Cy_f * tan_af * (vx_f / (norm_1(v_wf) + 1e-3));
        SX sigma_xf = Cx_f * (v_wf - vx_f) / (norm_1(v_wf) + 1e-3);
        SX gamma_f = norm_2(SX::vertcat({sigma_xf + reg_sign(sigma_xf) * 10, sigma_yf + reg_sign(sigma_yf) * 10}));

        SX gamma_sf = 3 * mu * FZf;
        SX f2f = mu * FZf;
        SX f1f = gamma_f - pow(gamma_f, 2) / gamma_sf + (1 / 3.0) * pow(gamma_f, 3) / pow(gamma_sf, 2);

        SX Fs_f = (1 - heaviside(gamma_f - gamma_sf + 1e-5)) * f1f + heaviside(gamma_f - gamma_sf + 1e-5) * f2f;
        FXf = (sigma_xf / (gamma_f + reg_sign(gamma_f) * 1e-3)) * Fs_f;
        FYf = (sigma_yf / (gamma_f + reg_sign(gamma_f) * 1e-3)) * Fs_f;
    }
    else if(tire == LINEAR)
    {
        /** longitudinal slip */
        SX k_f = (v_wf - vx_f) / (vx_f);     // should regularize here

        FXf = Cx_f * k_f;
        FYf = -Cy_f * atan2(vy_f, abs(vx_f));
    }

    SX FORCES = SX::vertcat({FXr, FYr, FXf, FYf});
    //SX FORCES = SX::vertcat({sigma_xf, sigma_yf, gamma_f, f1f, FXf, FYf});
    //SX sym_jac = SX::jacobian(FORCES, FXr_req);
    TraceFunction = Function("trace",{state, control},{FORCES});

    /** Drag force */
    /** compensate rollresist in the controller */
    Fdrag = sign(vx) * (rollResist + rollResistSpeed * vx + Cxx * pow(vx, 2));

    /** dynamic equations */
    SX vx_dot    = omega * vy + (FXr + FXf * cos(phi) - FYf * sin(phi) - Fdrag) / m;
    SX vy_dot    = -omega * vx + (FYf * cos(phi) + FXf * sin(phi) + FYr) / m;
    SX omega_dot = (-FYr * b + (FYf * cos(phi) + FXf * sin(phi)) * a) / Iz;

    SX x_dot     = vx * cos(theta) - vy * sin(theta);
    SX y_dot     = vx * sin(theta) + vy * cos(theta);
    SX theta_dot = omega;
    SX owf_dot   = r_f * (FXf_req - FXf) / fwIz;
    SX owr_dot   = r_r * (FXr_req - FXr) / rwIz;

    Dynamics = SX::vertcat({vx_dot, vy_dot, omega_dot, x_dot, y_dot, theta_dot, owf_dot, owr_dot});
    /** Dynamic equations */
    NumDynamics = Function("Dynamics", {state, control}, {Dynamics});

    /** define output mapping */
    SX H = SX::zeros(3,8);
    H(0,3) = 1; H(1,4) = 1; H(2,5) = 1;
    OutputMap = Function("Map",{state}, {SX::vertcat({x,y,theta})});
}

MobileRobot::MobileRobot()
{
    /** model constants [DevBot2]*/
    double g    = 9.81;    // gravity acceleration

    /** car parameters */
    MobileRobotProperties props; // default car properties
    double m    = props.mass;
    double b    = props.b;
    double a    = props.a;
    double Iz   = props.Iz;
    double r_f  = props.r_f;
    double r_r  = props.r_r;
    double fwIz = props.fwIz;
    double rwIz = props.rwIz;
    double Cy_r = props.Cy_r;
    double Cx_r = props.Cx_r;
    double Cy_f = props.Cy_f;
    double Cx_f = props.Cx_f;
    double mu   = props.mu;
    double rollResist = props.rollResist;
    double rollResistSpeed = props.rollResistSpeed;
    double Cxx = props.Cxx;
    TIRE_MODEL tire = LINEAR;


    /** state definition */
    SX vx    = SX::sym("vx");
    SX vy    = SX::sym("vy");
    SX omega = SX::sym("omega");

    SX x     = SX::sym("x");
    SX y     = SX::sym("y");
    SX theta = SX::sym("theta");

    SX owf   = SX::sym("owf");
    SX owr   = SX::sym("owr");
    state = SX::vertcat({vx, vy, omega, x, y, theta, owf, owr});

    /** controls definition */
    SX phi       = SX::sym("phi");
    SX FXr_req   = SX::sym("Fxr_req");
    SX FXf_req   = SX::sym("FXf_req");
    SX control   = SX::vertcat({phi, FXf_req, FXr_req});

    /** tire forces */
    SX FXf, FXr, FYf, FYr, Fdrag;

    /** vertical force acting on the wheels */
    SX FZr, FZf;

    /** REAR wheel force */
    SX v_wr = r_r * owr + 0.01;
    SX vx_r = vx + 0.01;
    SX vy_r = vy - omega * b;
    FZr = m * g * a / (a + b);

    if(tire == BRUSH)
    {
        SX tan_ar = tan( vy_r / (vx + reg_sign(vx) * 1e-3));   //regularization
        SX sigma_yr = -Cy_r * tan_ar * (vx / (norm_1(v_wr) + 1e-3));   // "regularization"
        SX sigma_xr = Cx_r * (v_wr - vx_r) / (norm_1(v_wr) + 1e-3);     //question Vsx sign
        SX gamma_r = norm_2(SX::vertcat({sigma_xr + reg_sign(sigma_xr) * 10, sigma_yr + reg_sign(sigma_yr) * 10}));

        /** Brush model of the tire */
        SX gamma_s = 3 * mu * FZr;
        SX f2 = mu * FZr;
        SX f1 = gamma_r - pow(gamma_r,2) / gamma_s + (1 / 3.0) * pow(gamma_r, 3) / pow(gamma_s, 2);

        SX Fs_r = (1 - heaviside(gamma_r - gamma_s + 1e-5)) * f1 + heaviside(gamma_r - gamma_s + 1e-5) * f2;
        FXr = (sigma_xr / (gamma_r)) * Fs_r;
        FYr = (sigma_yr / (gamma_r)) * Fs_r;
    }
    else if(tire == LINEAR)
    {
        /** longitudinal slip */
        SX k_r = (v_wr - vx_r) / (vx_r + 0.01);     // should regularize here

        FXr = Cx_r * k_r;
        FYr = -Cy_r * atan2(vy_r, abs(vx_r));                    //atan(tan_ar);
    }

    /** FRONT wheel force */
    SX v_wf = r_f * owf + 0.01;                                         // "regularization"
    SX vx_f = cos(-phi) * vx - sin(-phi) * (vy + omega * a) + 0.01;
    SX vy_f = sin(-phi) * vx + cos(-phi) * (vy + omega * a);
    FZf = m * g * b / (a + b);


    if(tire == BRUSH)
    {
        SX tan_af = tan( vy_f / (vx_f + reg_sign(vx_f) * 1e-3));
        SX sigma_yf = -Cy_f * tan_af * (vx_f / (norm_1(v_wf) + 1e-3));
        SX sigma_xf = Cx_f * (v_wf - vx_f) / (norm_1(v_wf) + 1e-3);
        SX gamma_f = norm_2(SX::vertcat({sigma_xf + reg_sign(sigma_xf) * 10, sigma_yf + reg_sign(sigma_yf) * 10}));

        SX gamma_sf = 3 * mu * FZf;
        SX f2f = mu * FZf;
        SX f1f = gamma_f - pow(gamma_f, 2) / gamma_sf + (1 / 3.0) * pow(gamma_f, 3) / pow(gamma_sf, 2);

        SX Fs_f = (1 - heaviside(gamma_f - gamma_sf + 1e-5)) * f1f + heaviside(gamma_f - gamma_sf + 1e-5) * f2f;
        FXf = (sigma_xf / (gamma_f + reg_sign(gamma_f) * 1e-3)) * Fs_f;
        FYf = (sigma_yf / (gamma_f + reg_sign(gamma_f) * 1e-3)) * Fs_f;
    }
    else if(tire == LINEAR)
    {
        /** longitudinal slip */
        SX k_f = (v_wf - vx_f) / (vx_f);     // should regularize here

        FXf = Cx_f * k_f;
        FYf = -Cy_f * atan2(vy_f, abs(vx_f));
    }

    SX FORCES = SX::vertcat({FXr, FYr, FXf, FYf});
    //SX FORCES = SX::vertcat({sigma_xf, sigma_yf, gamma_f, f1f, FXf, FYf});
    //SX sym_jac = SX::jacobian(FORCES, FXr_req);
    TraceFunction = Function("trace",{state, control},{FORCES});

    /** Drag force */
    /** compensate rollresist in the controller */
    Fdrag = sign(vx) * (rollResist + rollResistSpeed * vx + Cxx * pow(vx, 2));

    /** dynamic equations */
    SX vx_dot    = omega * vy + (FXr + FXf * cos(phi) - FYf * sin(phi) - Fdrag) / m;
    SX vy_dot    = -omega * vx + (FYf * cos(phi) + FXf * sin(phi) + FYr) / m;
    SX omega_dot = (-FYr * b + (FYf * cos(phi) + FXf * sin(phi)) * a) / Iz;

    SX x_dot     = vx * cos(theta) - vy * sin(theta);
    SX y_dot     = vx * sin(theta) + vy * cos(theta);
    SX theta_dot = omega;
    SX owf_dot   = r_f * (FXf_req - FXf) / fwIz;
    SX owr_dot   = r_r * (FXr_req - FXr) / rwIz;

    Dynamics = SX::vertcat({vx_dot, vy_dot, omega_dot, x_dot, y_dot, theta_dot, owf_dot, owr_dot});
    /** Dynamic equations */
    NumDynamics = Function("Dynamics", {state, control}, {Dynamics});

    /** define output mapping */
    SX H = SX::zeros(3,8);
    H(0,3) = 1; H(1,4) = 1; H(2,5) = 1;
    OutputMap = Function("Map",{state}, {SX::vertcat({x,y,theta})});
}
