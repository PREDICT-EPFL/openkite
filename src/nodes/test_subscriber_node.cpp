#include "ros/ros.h"
#include "std_msgs/String.h"
#include "eigen3/Eigen/Dense"


void chatterCallback(const std_msgs::String::ConstPtr& msg)
{
  ROS_INFO("PONG to : [%s]", msg->data.c_str());
}

void controlCb(const std_msgs::String::ConstPtr& msg)
{
    ROS_INFO("CONTROL RECEIVED : [%s]", msg->data.c_str());
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "test_subscriber");
  ros::NodeHandle n;

  ros::Subscriber sub = n.subscribe("chatter", 1000, chatterCallback);
  ros::Subscriber sub_ctl = n.subscribe("control", 10, controlCb);

  ros::spin();

  return 0;
}
