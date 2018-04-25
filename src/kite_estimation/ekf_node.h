#ifndef EKF_NODE_H
#define EKF_NODE_H

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "geometry_msgs/PoseStamped.h"
#include "sensor_msgs/MultiDOFJointState.h"
#include "std_msgs/Int16MultiArray.h"

#include "boost/thread/mutex.hpp"
#include "boost/circular_buffer.hpp"

#include "kiteEKF.h"
#include "sstream"

class KiteEKF_Node
{
public:
    KiteEKF_Node(const ros::NodeHandle &_nh, const casadi::Function &Integrator,
                                             const casadi::Function &Jacobian );
    virtual ~KiteEKF_Node(){}

    sensor_msgs::MultiDOFJointState kite_state;
    boost::circular_buffer<geometry_msgs::PoseStamped> measurements;
    ros::Time last_measurement_arrived;

    void filterCallback(const geometry_msgs::PoseStamped::ConstPtr &msg);
    void controlCallback(const std_msgs::Int16MultiArray::ConstPtr &msg);

    void estimate(geometry_msgs::PoseStamped _pose);
    void estimate();
    void publish();

    void initialize();
    bool initialized(){return m_initialized;}

    void set_brf_offset(const casadi::DM &offset){brf_offset = offset;}
    void set_brf_rotation(const casadi::DM &rotation){brf_rotation = rotation;}

    boost::mutex m_mutex;

private:
    casadi::DM control;
    casadi::DM brf_offset;
    casadi::DM brf_rotation;

    ros::Publisher state_pub;
    ros::Subscriber control_sub;
    ros::Subscriber pose_sub;

    /** handle instance to access node params */
    std::shared_ptr<ros::NodeHandle> nh;
    /** EKF instance */
    std::shared_ptr<KiteEKF> filter;
    bool m_initialized;

    casadi::DM convertToDM(const geometry_msgs::PoseStamped &_value);
    casadi::DMVector convertToDM(const sensor_msgs::MultiDOFJointState &_value);
    casadi::DM optitrack2world(const casadi::DM &opt_pose);
};


#endif // EKF_NODE_H
