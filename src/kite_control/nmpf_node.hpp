#ifndef NMPF_NODE_HPP
#define NMPF_NODE_HPP

#include "ros/ros.h"
#include "sensor_msgs/MultiDOFJointState.h"
#include "std_msgs/Int32MultiArray.h"
#include "openkite/aircraft_controls.h"
#include "geometry_msgs/PoseStamped.h"
#include "openkite/mpc_diagnostic.h"

#include "boost/thread/mutex.hpp"
#include "kiteNMPF.h"
#include "integrator.h"
#include "nmpf.hpp"
#include "mobile_robot.h"

using namespace casadi;

/** polynomial helper function */
SX polynomial(const SX &coef, const SX &x0, const SX &x)
{
    SX poly = SX::vertcat(SXVector{pow(x-x0, 3), pow(x-x0, 2), x-x0, 1});
    return SX::dot(coef, poly);
}

struct Path
{
    SXVector operator()(const SXVector &arg)
    {
        SX breaks = DM({0, 1.4225, 2.8451, 4.2676, 5.6902, 7.1127, 8.5353, 9.9578, 11.3803});
        SX coeffs_x = DM::horzcat(DMVector{DM({-0.0049,    0.0510,   -0.0261,   -5.0000}),
                                         DM({-0.0049,    0.0300,    0.0890,   -4.9482}),
                                         DM({0.0555,    0.0090,    0.1445,   -4.7750}),
                                         DM({-0.0635,    0.2459,    0.5071,   -4.3915}),
                                         DM({-0.0173,   -0.0252,    0.8211,   -3.3554}),
                                         DM({0.0212,   -0.0989,    0.6446,   -2.2881}),
                                         DM({-0.0003,   -0.0082,    0.4923,   -1.5100}),
                                         DM({-0.0003,   -0.0093,    0.4674,   -0.8270}) });

        SX coeffs_y = DM::horzcat(DMVector{DM({0.0005,   -0.0047,    1.0048,   -0.0000}),
                                           DM({0.0005,   -0.0025,    0.9945,    1.4213}),
                                           DM({-0.0148,   -0.0002,    0.9908,    2.8326}),
                                           DM({-0.1311,   -0.0632,    0.9007,    4.1992}),
                                           DM({0.1598,   -0.6226,   -0.0749,    4.9753}),
                                           DM({-0.0221,    0.0594,   -0.8761,    4.0689}),
                                           DM({0.0071,   -0.0351,   -0.8415,    2.8790}),
                                           DM({0.0071,   -0.0050,   -0.8984,    1.6313}) });

        SX x = SX::sym("x");
        coeffs_x = coeffs_x.T();
        coeffs_y = coeffs_y.T();

        SX path_x = 0;
        SX path_y = 0;
        for(int i = 0; i <= 7; ++i)
        {
            SX pix = polynomial(coeffs_x(i, Slice(0, 4)).T(), SX(breaks[i]), x);
            path_x += pix * (heaviside(x - breaks[i] + 1e-5) - heaviside(x - breaks[i+1])); // 1e-5 - Heaviside correction

            SX piy = polynomial(coeffs_y(i, Slice(0, 4)).T(), SX(breaks[i]), x);
            path_y += piy * (heaviside(x - breaks[i] + 1e-5) - heaviside(x - breaks[i+1])); // 1e-5 - Heaviside correction
        }
        SX Path = SX::vertcat({path_x, path_y});
        Function path = Function("Path", {x}, {Path});
        return path(arg);
    }
};

using NMPF = polympc::nmpf<MobileRobot, Path, 3, 2>;

class KiteNMPF_Node
{
public:
    KiteNMPF_Node(const ros::NodeHandle &_nh);
    virtual ~KiteNMPF_Node(){}

    ros::Time last_computed_control;
    void filterCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg);

    //void compute_control(const geometry_msgs::PoseStamped &_pose);
    void compute_control();
    void publish();
    void publish_trajectory();
    void publish_mpc_diagnostic();

    void initialize(){m_initialized = true;}
    bool is_initialized(){return m_initialized;}

    boost::mutex m_mutex;
    double comp_time_ms;

private:
    casadi::DM control;
    casadi::DM kite_state;

    ros::Publisher  control_pub;
    ros::Publisher  traj_pub;
    ros::Publisher  diagnostic_pub;
    ros::Subscriber state_sub;

    /** handle instance to access node params */
    std::shared_ptr<ros::NodeHandle> nh;
    /** NMPF instance */
    std::shared_ptr<NMPF> controller;
    std::shared_ptr<ODESolver> solver;

    bool m_initialized;
    double transport_delay;

    casadi::DM convertToDM(const sensor_msgs::MultiDOFJointState &_value);
};


#endif // NMPF_NODE_HPP

