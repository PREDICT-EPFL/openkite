#ifndef CONTROL_PROXY_NODE_HPP
#define CONTROL_PROXY_NODE_HPP

#include "ros/ros.h"
#include "std_msgs/Int32MultiArray.h"
#include "std_msgs/Int16MultiArray.h"
#include "openkite/aircraft_controls.h"
#include "boost/thread/mutex.hpp"

class ControlProxyNode
{
public:
    ControlProxyNode(const ros::NodeHandle &_nh);
    virtual ~ControlProxyNode(){}

    void publish(){ proxy_pub.publish(servo_msg); }
    void set_servos(const int &thrust, const int &elevator,
                    const int &rudder, const int &ailerons)
    {
        servo_msg.data.clear();
        servo_msg.data.push_back( static_cast<uint16_t>(thrust) );
        servo_msg.data.push_back( static_cast<uint16_t>(rudder) );
        servo_msg.data.push_back( static_cast<uint16_t>(elevator) );
        servo_msg.data.push_back( static_cast<uint16_t>(ailerons) );
    }

private:
    std_msgs::Int16MultiArray servo_msg;
    void controlCallback(const openkite::aircraft_controls::ConstPtr& msg);
    ros::Subscriber proxy_sub;
    ros::Publisher proxy_pub;
    std::shared_ptr<ros::NodeHandle> nh;

    boost::mutex m_mutex;
};

#endif // CONTROL_PROXY_NODE_HPP

