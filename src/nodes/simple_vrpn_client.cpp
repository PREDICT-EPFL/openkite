#include "simple_vrpn_client.h"


#include "tf2/LinearMath/Quaternion.h"
#include "tf2/LinearMath/Matrix3x3.h"
#include "tf2_ros/transform_broadcaster.h"

#include <vector>
#include <unordered_set>

namespace
{
  std::unordered_set<std::string> name_blacklist_({"VRPN Control"});
}

namespace vrpn_client_ros
{

  VrpnTrackerRos::VrpnTrackerRos(std::string tracker_name, ConnectionPtr connection, ros::NodeHandle nh)
  {
    tracker_remote_ = std::make_shared<vrpn_Tracker_Remote>(tracker_name.c_str(), connection.get());
    init(tracker_name, nh, false);
  }

  VrpnTrackerRos::VrpnTrackerRos(std::string tracker_name, std::string host, ros::NodeHandle nh)
  {
    std::string tracker_address;
    tracker_address = tracker_name + "@" + host;
    ROS_INFO("Tracker address: %s \n", tracker_address.c_str());
    tracker_remote_ = std::make_shared<vrpn_Tracker_Remote>(tracker_address.c_str());
    init(tracker_name, nh, true);
  }

  void VrpnTrackerRos::init(std::string tracker_name, ros::NodeHandle nh, bool create_mainloop_timer)
  {
    ROS_INFO_STREAM("Creating new tracker " << tracker_name);

    tracker_remote_->register_change_handler(this, &VrpnTrackerRos::handle_pose);
    //tracker_remote_->shutup = true;

    std::string error;
    if (!ros::names::validate(tracker_name, error))
    {
      ROS_ERROR_STREAM("Invalid tracker name " << tracker_name << ", not creating topics : " << error);
      return;
    }

    this->tracker_name = tracker_name;

    output_nh_ = ros::NodeHandle(nh, tracker_name);

    std::string frame_id;
    nh.param<std::string>("frame_id", frame_id, "world");
    nh.param<bool>("use_server_time", use_server_time_, false);
    nh.param<bool>("broadcast_tf", broadcast_tf_, false);
    nh.param<bool>("process_sensor_id", process_sensor_id_, false);

    pose_msg_.header.frame_id = frame_id;

    if (create_mainloop_timer)
    {
      double update_frequency;
      nh.param<double>("update_frequency", update_frequency, 100.0);
      mainloop_timer = nh.createTimer(ros::Duration(1 / update_frequency),
                                      boost::bind(&VrpnTrackerRos::mainloop, this));
    }
  }

  VrpnTrackerRos::~VrpnTrackerRos()
  {
    ROS_INFO_STREAM("Destroying tracker " << pose_msg_.header.frame_id);
    tracker_remote_->unregister_change_handler(this, &VrpnTrackerRos::handle_pose);
  }

  void VrpnTrackerRos::mainloop()
  {
      tracker_remote_->mainloop();
  }

  void VRPN_CALLBACK VrpnTrackerRos::handle_pose(void *userData, const vrpn_TRACKERCB tracker_pose)
  {

    VrpnTrackerRos *tracker = static_cast<VrpnTrackerRos *>(userData);

    ros::Publisher *pose_pub;
    //std::size_t sensor_index(0);
    ros::NodeHandle nh = tracker->output_nh_;

    if (tracker->process_sensor_id_)
    {
      //sensor_index = static_cast<std::size_t>(tracker_pose.sensor);
      nh = ros::NodeHandle(tracker->output_nh_, std::to_string(tracker_pose.sensor));
    }

    pose_pub = &(tracker->pose_pubs_);

    if (pose_pub->getTopic().empty())
    {
      *pose_pub = nh.advertise<geometry_msgs::PoseStamped>("pose", 1);
    }

    //ROS_INFO("POS [%f %f %f] ATT [%f %f %f %f] \n", tracker_pose.pos[0], tracker_pose.pos[1], tracker_pose.pos[2],
    //        tracker_pose.quat[0], tracker_pose.quat[1], tracker_pose.quat[2], tracker_pose.quat[3]);

    if (pose_pub->getNumSubscribers() > 0)
    {
      if (tracker->use_server_time_)
      {
        tracker->pose_msg_.header.stamp.sec = tracker_pose.msg_time.tv_sec;
        tracker->pose_msg_.header.stamp.nsec = tracker_pose.msg_time.tv_usec * 1000;
      }
      else
      {
        tracker->pose_msg_.header.stamp = ros::Time::now();
      }

      tracker->pose_msg_.pose.position.x = tracker_pose.pos[0];
      tracker->pose_msg_.pose.position.y = tracker_pose.pos[1];
      tracker->pose_msg_.pose.position.z = tracker_pose.pos[2];

      tracker->pose_msg_.pose.orientation.x = tracker_pose.quat[0];
      tracker->pose_msg_.pose.orientation.y = tracker_pose.quat[1];
      tracker->pose_msg_.pose.orientation.z = tracker_pose.quat[2];
      tracker->pose_msg_.pose.orientation.w = tracker_pose.quat[3];

      //ROS_INFO_STREAM("optitrack: caught a measurement");

      pose_pub->publish(tracker->pose_msg_);
    }
  }

  //-----------------------------------------------------------------//

  VrpnClientRos::VrpnClientRos(ros::NodeHandle nh, ros::NodeHandle private_nh)
  {
    output_nh_ = private_nh;

    host_ = getHostStringFromParams(private_nh);

    ROS_INFO_STREAM("Connecting to VRPN server at " << host_);
    connection_ = std::shared_ptr<vrpn_Connection>(vrpn_get_connection_by_name(host_.c_str()));
    ROS_INFO("Connection established");

    double update_frequency;
    private_nh.param<double>("update_frequency", update_frequency, 50.0);
    mainloop_timer = nh.createTimer(ros::Duration(1 / update_frequency), boost::bind(&VrpnClientRos::mainloop, this));

    double refresh_tracker_frequency;
    private_nh.param<double>("refresh_tracker_frequency", refresh_tracker_frequency, 0.0);

    // use only one tracker
    std::string param_tracker_name_;
    //private_nh.getParam("trackers", param_tracker_name_);
    private_nh.param<std::string>("trackers",param_tracker_name_,"Kite");
    trackers_.insert(std::make_pair(param_tracker_name_,
                                    std::make_shared<VrpnTrackerRos>(param_tracker_name_, connection_, output_nh_)));
  }

  std::string VrpnClientRos::getHostStringFromParams(ros::NodeHandle host_nh)
  {
    std::stringstream host_stream;
    std::string server;
    int port;

    host_nh.param<std::string>("server", server, "192.168.0.249");
    host_nh.param<int>("port", port, DEFAULT_OPTITRACK_PORT);
    host_stream << "tcp://" << server << ":" << port;

    ROS_INFO("Connecting to: %s \n", host_stream.str().c_str());

    return host_stream.str();
  }

  void VrpnClientRos::mainloop()
  {
    connection_->mainloop();

    if (!connection_->doing_okay())
    {
      ROS_WARN("VRPN connection is not 'doing okay'");
    }
    for (TrackerMap::iterator it = trackers_.begin(); it != trackers_.end(); ++it)
    {
        it->second->mainloop();
    }
  }

}  // namespace vrpn_client_ros
