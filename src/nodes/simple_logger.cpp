#include "ros/ros.h"
#include "ros/time.h"
#include "std_msgs/Int32MultiArray.h"
#include "geometry_msgs/PoseStamped.h"
#include "sensor_msgs/MultiDOFJointState.h"

#include "iostream"
#include "fstream"
#include "string.h"


class SimpleLogger
{
public:
    SimpleLogger();
    virtual ~SimpleLogger()
    {
        /** close loggers */
        control_log.close();
        pose_log.close();
        state_log.close();
    }

    void controlCallback(const std_msgs::Int32MultiArray::ConstPtr& msg);
    void poseCallback(const geometry_msgs::PoseStamped::ConstPtr& msg);
    void stateCallback(const sensor_msgs::MultiDOFJointState::ConstPtr& msg);

private:
    std::ofstream control_log;
    std::ofstream pose_log;
    std::ofstream state_log;
};

SimpleLogger::SimpleLogger()
{
    /** create logger files */
    control_log.open("kite_control.log");
    pose_log.open("kite_pose.log");
    state_log.open("kite_state.log");

    if(control_log.is_open())
        std::cout << "Kite control log file created in the current directory \n";

    if(pose_log.is_open())
        std::cout << "Kite pose log file created in the current directory \n";

    if(state_log.is_open())
        std::cout << "Kite state estimation log file created in the current directory \n";
}

void SimpleLogger::controlCallback(const std_msgs::Int32MultiArray::ConstPtr& msg)
{
    /** @badcode: it's not recommended to precess data in callbacks */
    control_log << std::fixed << std::setprecision(8) << ros::Time::now().toSec() << " ";
    for (unsigned i = 0; i < 4; ++i)
    {
        control_log << msg->data[i] << " ";
    }
    control_log << "\n";
}

void SimpleLogger::poseCallback(const geometry_msgs::PoseStamped::ConstPtr &msg)
{
    pose_log << std::fixed << std::setprecision(8) << msg->header.stamp << " ";
    pose_log << msg->pose.position.x << " " << msg->pose.position.y << " " << msg->pose.position.z << " "
             << msg->pose.orientation.w << " " << msg->pose.orientation.x << " " << msg->pose.orientation.y << " " << msg->pose.orientation.z << "\n";
}

/** state estimation logger callback */
void SimpleLogger::stateCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg)
{
    if (msg->transforms.empty() || msg->twist.empty())
    {
        std::cout << "Message came but empty \n";
    }
    else
    {
        state_log << std::fixed << std::setprecision(8) << msg->header.stamp << " ";
        state_log << msg->twist[0].linear.x << " " << msg->twist[0].linear.y << " " << msg->twist[0].linear.z << " "
                                           << msg->twist[0].angular.x << " " << msg->twist[0].angular.y << " " << msg->twist[0].angular.z << " "
                                           << msg->transforms[0].translation.x << " " << msg->transforms[0].translation.y << " " << msg->transforms[0].translation.z << " "
                                           << msg->transforms[0].rotation.w << " " << msg->transforms[0].rotation.x
                                           << " " << msg->transforms[0].rotation.y << " " << msg->transforms[0].rotation.z << "\n";
    }
}

int main(int argc, char **argv)
{
    SimpleLogger logger;

    ros::init(argc, argv, "simple_logger");
    ros::NodeHandle n;

    ros::Subscriber sub_ctrl = n.subscribe("chatter", 100, &SimpleLogger::controlCallback, &logger);
    ros::Subscriber sub_pose = n.subscribe("optitrack_client/RigidBody1/pose", 100, &SimpleLogger::poseCallback, &logger);
    ros::Subscriber sub_state = n.subscribe("kite_state", 100, &SimpleLogger::stateCallback, &logger);

    ros::spin();
    return 0;
}
