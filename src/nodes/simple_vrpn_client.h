#ifndef SIMPLE_VRPN_CLIENT_H
#define SIMPLE_VRPN_CLIENT_H

#include "ros/ros.h"
#include "geometry_msgs/PoseStamped.h"

#include <vrpn_Tracker.h>
#include <vrpn_Connection.h>
#include <string>

#define DEFAULT_OPTITRACK_PORT 3883


namespace vrpn_client_ros
{

  typedef std::shared_ptr<vrpn_Connection> ConnectionPtr;
  typedef std::shared_ptr<vrpn_Tracker_Remote> TrackerRemotePtr;

  class VrpnTrackerRos
  {
  public:

    typedef std::shared_ptr<VrpnTrackerRos> Ptr;
    /**
     * Create and initialize VrpnTrackerRos using an existing underlying VRPN connection object. The underlying
     * connection object is responsible for calling the tracker's mainloop.
     */
    VrpnTrackerRos(std::string tracker_name, ConnectionPtr connection, ros::NodeHandle nh);

    /**
     * Create and initialize VrpnTrackerRos, creating a new connection to tracker_name@host. This constructor will
     * register timer callbacks on nh to call mainloop.
     */
    VrpnTrackerRos(std::string tracker_name, std::string host, ros::NodeHandle nh);

    ~VrpnTrackerRos();

    /**
     * Call mainloop of underlying vrpn_Tracker_Remote
     */
    void mainloop();

  private:
    TrackerRemotePtr tracker_remote_;
    ros::Publisher pose_pubs_;
    ros::NodeHandle output_nh_;
    bool use_server_time_, broadcast_tf_, process_sensor_id_;
    std::string tracker_name;

    ros::Timer mainloop_timer;

    geometry_msgs::PoseStamped pose_msg_;

    void init(std::string tracker_name, ros::NodeHandle nh, bool create_mainloop_timer);

    static void VRPN_CALLBACK handle_pose(void *userData, const vrpn_TRACKERCB tracker_pose);
  };

  class VrpnClientRos
  {
  public:

    typedef std::shared_ptr<VrpnClientRos> Ptr;
    typedef std::map<std::string, VrpnTrackerRos::Ptr> TrackerMap;

    /**
     * Create and initialize VrpnClientRos object in the private_nh namespace.
     */
    VrpnClientRos(ros::NodeHandle nh, ros::NodeHandle private_nh);

    static std::string getHostStringFromParams(ros::NodeHandle host_nh);

    /**
     * Call mainloop of underlying VRPN connection and all registered VrpnTrackerRemote objects.
     */
    void mainloop();

    /**
     * Examine vrpn_Connection's senders and create new trackers as necessary.
     */
    void updateTrackers();

  private:
    std::string host_;
    ros::NodeHandle output_nh_;

    /**
     * Underlying VRPN connection object
     */
    ConnectionPtr connection_;

    /**
     * Map of registered trackers, accessible by name
     */
    TrackerMap trackers_;

    ros::Timer refresh_tracker_timer_, mainloop_timer;
  };
}  // namespace vrpn_client_ros

#endif // SIMPLE_VRPN_CLIENT_H
