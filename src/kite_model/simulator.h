#ifndef SIMULATOR_H
#define SIMULATOR_H

//#include "stochastic_kite.hpp"
#include "mobile_robot.h"
#include "kinematic_robot.hpp"
#include "kite.h"
#include "integrator.h"
#include "ros/ros.h"
#include "openkite/aircraft_controls.h"
#include "sensor_msgs/MultiDOFJointState.h"
#include "geometry_msgs/PoseStamped.h"


class Simulator
{
public:
    Simulator(const ODESolver &object, const ros::NodeHandle &nh);
    virtual ~Simulator(){}
    void simulate();

    casadi::DM getState(){return state;}
    void setControls(const casadi::DM &_control){controls = _control;}
    //casadi::DM getPose(){return state(casadi::Slice(6,13));}

    void publish_state();

    bool is_initialized(){return initialized;}
    void initialize(const casadi::DM &_init_value){state = _init_value; initialized = true;}

private:
    std::shared_ptr<ODESolver> m_object;
    std::shared_ptr<ros::NodeHandle> m_nh;

    ros::Subscriber control_sub;
    ros::Publisher  state_pub;
    ros::Publisher  pose_pub;

    casadi::DM      controls;
    casadi::DM      state;

    void controlCallback(const openkite::aircraft_controls::ConstPtr &msg);
    double sim_rate;
    sensor_msgs::MultiDOFJointState msg_state;

    bool initialized;
};

#endif // SIMULATOR_H
