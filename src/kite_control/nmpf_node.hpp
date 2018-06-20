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

class KiteNMPF_Node
{
public:
    KiteNMPF_Node(const ros::NodeHandle &_nh, const SimpleKinematicKiteProperties &kite_props,
                                              const AlgorithmProperties &algo_props );
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
    std::shared_ptr<KiteNMPF> controller;
    std::shared_ptr<ODESolver> solver;

    bool m_initialized;
    double transport_delay;

    casadi::DM convertToDM(const sensor_msgs::MultiDOFJointState &_value);
};


#endif // NMPF_NODE_HPP
