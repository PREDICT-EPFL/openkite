#include "kinematic_robot.hpp"

using namespace casadi;

KinematicRobot::KinematicRobot(const MobileRobotProperties &props)
{

}

KinematicRobot::KinematicRobot()
{
    /** model constants [DevBot2]*/
    double g    = 9.81;    // gravity acceleration

    /** car parameters */
    MobileRobotProperties props; // default car properties
    double m    = props.mass;
    double b    = props.b;
    double a    = props.a;
    double Iz   = props.Iz;
    double Cy_r = props.Cy_r;
    double Cy_f = props.Cy_f;
    double rollResist = props.rollResist;
    double Cxx = props.Cxx;

    /** state and control description */
    SX vx     = SX::sym("vx");
    SX vy     = SX::sym("vy");
    SX omega = SX::sym("omega");
    SX x     = SX::sym("x");
    SX y     = SX::sym("y");
    SX theta = SX::sym("theta");
    state = SX::vertcat({vx, vy, omega, x, y, theta});

    SX phi = SX::sym("phi");
    SX FXr_req = SX::sym("FXr_req");
    SX FXf_req = SX::sym("FXf_req");
    control = SX::vertcat({phi, FXf_req, FXr_req});

    /** linear tire model */
    /** tire forces */
    SX FXf, FXr, FYf, FYr, Fdrag;

    /** vertical force acting on the wheels */
    SX FZr, FZf;

    /** REAR wheel force */
    SX vx_r = vx;
    SX vy_r = vy - omega * b;
    FZr = m * g * a / (a + b);
    //lateral slip force
    FYr = -Cy_r * atan2(vy_r, abs(vx_r)+0.01);

    /** FRONT wheel force */
    SX vx_f = cos(-phi) * vx - sin(-phi) * (vy + omega * a) ;
    SX vy_f = sin(-phi) * vx + cos(-phi) * (vy + omega * a);
    FZf = m * g * b / (a + b);
    //lateral slip force
    FYf = -Cy_f * atan2(vy_f, abs(vx_f)+0.01);

    FXr = FXr_req;
    FXf = FXf_req;

    SX FORCES = SX::vertcat({FXr, FYr, FXf, FYf});
    //SX FORCES = SX::vertcat({sigma_xf, sigma_yf, gamma_f, f1f, FXf, FYf});
    //SX sym_jac = SX::jacobian(FORCES, FXr_req);
    TraceFunction = Function("trace",{state, control},{FORCES});

    Fdrag = sign(vx) * (rollResist + Cxx * pow(vx, 2));

    /** dynamic equations */
    SX vx_dot    = omega * vy + (FXr_req + FXf_req * cos(phi) - FYf * sin(phi) ) / m;
    SX vy_dot    = -omega * vx + (FYf * cos(phi) + FXf_req * sin(phi) + FYr) / m;
    SX omega_dot = (-FYr * b + (FYf * cos(phi) + FXf_req * sin(phi)) * a) / Iz;

    SX x_dot     = vx * cos(theta) - vy * sin(theta);
    SX y_dot     = vx * sin(theta) + vy * cos(theta);
    SX theta_dot = omega;


    Dynamics = SX::vertcat({vx_dot, vy_dot, omega_dot, x_dot, y_dot, theta_dot});
    /** Dynamic equations */
    NumDynamics = Function("Dynamics", {state, control}, {Dynamics});

    /** define output mapping */
    OutputMap = Function("Map",{state}, {SX::vertcat({x,y,theta})});
}
