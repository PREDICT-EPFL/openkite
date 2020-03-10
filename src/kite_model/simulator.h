#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "kite.h"
#include "integrator.h"
#include "ros/ros.h"
#include "sensor_msgs/MultiDOFJointState.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/Vector3Stamped.h"
#include "openkite/aircraft_controls.h"
#include <sensor_msgs/Joy.h>

class Simulator
{
public:
    Simulator(const ODESolver &odeSolver, const ros::NodeHandle &nh);
    virtual ~Simulator(){}
    void simulate();

    casadi::DM getState(){return state;}
    casadi::DM getPose(){return state(casadi::Slice(6,13));}

    void publish_state();
    void publish_pose();

    bool is_initialized(){return initialized;}
    void initialize(const casadi::DM &_init_value){state = _init_value; initialized = true;}

    void setNumericSpecNongravForce(const casadi::Function &_NumericSpecNongravForce) {m_NumericSpecNongravForce = _NumericSpecNongravForce;}

private:
    std::shared_ptr<ros::NodeHandle> m_nh;
    std::shared_ptr<ODESolver> m_odeSolver;
    casadi::Function m_NumericSpecNongravForce;

    ros::Subscriber control_sub;
    ros::Publisher  state_pub;
    ros::Publisher  pose_pub;
    ros::Publisher  accel_pub;

    casadi::DM      controls;
    casadi::DM      state;
    std::vector<double> specNongravForce;

//    void controlCallback(const openkite::aircraft_controls::ConstPtr &msg);
    void controlCallback(const sensor_msgs::JoyConstPtr &msg);
    double sim_rate;
    sensor_msgs::MultiDOFJointState msg_state;
    //geometry_msgs::Vector3Stamped msg_accel;

    bool initialized;
};

#endif // SIMULATOR_H
