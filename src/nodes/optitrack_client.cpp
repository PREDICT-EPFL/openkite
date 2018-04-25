#include "simple_vrpn_client.h"

int main(int argc, char **argv)
{
  ros::init(argc, argv, "optitrack_client");
  ros::NodeHandle nh, private_nh("~");
  vrpn_client_ros::VrpnClientRos client(nh, private_nh);
  ros::spin();
  return 0;
}
