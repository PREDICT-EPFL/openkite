#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include "sensor_msgs/MultiDOFJointState.h"
#include "geometry_msgs/Vector3.h"
#include "geometry_msgs/Quaternion.h"
#include "geometry_msgs/Point.h"
#include "geometry_msgs/Transform.h"

#include "kitemath.h"

using namespace casadi;
using namespace kmath;

struct Scale
{
    float x;
    float y;
    float z;

    Scale(const float &_x=1, const float &_y=1, const float &_z=1) : x(_x), y(_y), z(_z) {}
    virtual ~Scale(){}
};

struct Color
{
    float r;
    float g;
    float b;
    float transparency;

    Color(const float &_r=0, const float &_g=0, const float &_b=1, const float &_t=0) : r(_r), g(_g), b(_b), transparency(_t){}
    virtual ~Color(){}
};

struct MarkerProperties
{
    int         type;
    int         id;
    int         action;
    std::string frame_id;
    std::string ns;
    std::string mesh_resource;
    Scale       scale;
    Color       color;

    geometry_msgs::Point      position;
    geometry_msgs::Quaternion orientation;

    MarkerProperties() : type(visualization_msgs::Marker::SPHERE),
                         id(0),
                         action(visualization_msgs::Marker::ADD),
                         frame_id("/kite"),
                         ns("awe"),
                         scale(Scale()),
                         color(Color()),
                         mesh_resource("")
    {
        position.x = position.y = position.z = 0;
        orientation.x = orientation.y = orientation.z = 0;
        orientation.w = 1;
    }

    void configureMarker(visualization_msgs::Marker &marker);
};

void MarkerProperties::configureMarker(visualization_msgs::Marker &marker)
{
    marker.type = type;
    marker.id   = id;
    marker.ns   = ns;
    marker.header.frame_id = frame_id;
    marker.header.stamp    = ros::Time::now();
    marker.action = action;
    marker.color.r = color.r; marker.color.g = color.g; marker.color.b = color.b;
    marker.color.a = color.transparency;
    marker.scale.x = scale.x; marker.scale.y = scale.y; marker.scale.z = scale.z;
    marker.pose.position = position;
    marker.pose.orientation = orientation;
    if(type == visualization_msgs::Marker::MESH_RESOURCE)
        marker.mesh_resource    = mesh_resource;
}

class KiteVisualizer
{
public:
    KiteVisualizer(const ros::NodeHandle &_nh);
    virtual ~KiteVisualizer(){}

    geometry_msgs::Point getTranslation(){geometry_msgs::Vector3 point = kite_state.transforms[0].translation;
                                          geometry_msgs::Point p; p.x = point.x; p.y = point.y; p.z = point.z;
                                          return p;}

    geometry_msgs::Quaternion getRotation(){return kite_state.transforms[0].rotation;}
    bool is_initialized(){return initialized;}
    uint32_t getNumPublishers(){return state_sub.getNumPublishers();}

    DM world2rviz(const casadi::DM &pose);
    void state2marker(const casadi::DM &kite_pose, visualization_msgs::Marker &_marker);
    visualization_msgs::Marker getPose(){return m_kite_marker;}

    visualization_msgs::MarkerArray getOptimalTrajectory();
    visualization_msgs::Marker      getVirtualTrajectory();

    visualization_msgs::Marker      getKiteMarker(){return m_kite_marker;}
    visualization_msgs::Marker      getTetherMarker();

    void setPath(const Function &path){PathFunction = path;}
    bool tether_active;

private:
    std::shared_ptr<ros::NodeHandle> nh;
    ros::Subscriber state_sub;
    ros::Subscriber opt_traj_sub;

    void filterCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg);
    void controlCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg);

    sensor_msgs::MultiDOFJointState    kite_state;
    visualization_msgs::Marker         m_kite_marker;

    /** optimal trajectory visualisation */
    visualization_msgs::MarkerArray    kite_traj_markers;
    visualization_msgs::Marker         virtual_state_markers;
    visualization_msgs::Marker         tether_marker;

    DM   optimal_trajectory;
    Function PathFunction;
    bool initialized;

    DM convertToDM(const sensor_msgs::MultiDOFJointState &_value, const bool &with_virtual = false);
};

KiteVisualizer::KiteVisualizer(const ros::NodeHandle &_nh)
{
    nh = std::make_shared<ros::NodeHandle>(_nh);

    kite_state.twist.resize(1);
    kite_state.transforms.resize(1);
    kite_state.joint_names.resize(1);
    kite_state.header.frame_id = "kite";
    kite_state.joint_names[0] = "lox";

    std::string state_topic = "/kite_state";
    state_sub = nh->subscribe(state_topic, 100, &KiteVisualizer::filterCallback, this);
    opt_traj_sub = nh->subscribe("/opt_traj", 100, &KiteVisualizer::controlCallback, this);

    /** configure kite_marker */
    MarkerProperties m_props;
    m_props.id = 101;
    m_props.type = visualization_msgs::Marker::SPHERE;
    //m_props.type = visualization_msgs::Marker::MESH_RESOURCE;
    //m_props.mesh_resource = "package://openkite/meshes/kite_rviz.dae";
    m_props.scale = Scale(0.2, 0.2, 0.2);
    m_props.color = Color(1.0, 0.0, 0.0, 1.0);
    m_props.configureMarker(m_kite_marker);
    m_kite_marker.lifetime = ros::Duration();

    /** configure virtual state trajectory */
    m_props.id = 201;
    m_props.type  = visualization_msgs::Marker::SPHERE_LIST;
    m_props.color = Color(0.0, 0.0, 1.0, 0.5);
    m_props.configureMarker(virtual_state_markers);
    virtual_state_markers.lifetime = ros::Duration();

    /** configure tether marker */
    m_props.id = 301;
    m_props.type = visualization_msgs::Marker::LINE_STRIP;
    m_props.color = Color(1.0, 0.5, 0.0, 1.0);
    m_props.scale.x = 0.025;
    m_props.configureMarker(tether_marker);
    tether_marker.lifetime = ros::Duration();
    tether_marker.points.clear();
    geometry_msgs::Point p;
    p.x = 0.0;
    p.y = 0.0;
    p.z = 0.0;
    tether_marker.points.resize(2);
    tether_marker.points[0] = p;
    tether_marker.points[1] = p;

    initialized = false;
    tether_active = false;
}

void KiteVisualizer::filterCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg)
{
    /** transform kite pose to Rviz reference frame */
    /** @badcode: */
    DM pose = convertToDM(*msg);
    DM norm = DM::norm_2(pose(Slice(0,3)));
    if(norm.nonzeros()[0] >= 2.74)
        tether_active = true;
    else
        tether_active = false;

    //DM kite_state_dm = world2rviz(pose);
    state2marker(pose, m_kite_marker);
    initialized = true;
}

void KiteVisualizer::controlCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg)
{
    optimal_trajectory = convertToDM(*msg, true);
}

DM KiteVisualizer::convertToDM(const sensor_msgs::MultiDOFJointState &_value, const bool &with_virtual)
{
    uint num_poses = _value.transforms.size();
    int  size      = with_virtual ? 10 : 7;
    DM poses = DM::zeros(size, num_poses);

    for(uint i = 0; i < num_poses; ++i)
    {
        poses(0, i) = _value.transforms[i].translation.x;
        poses(1, i) = _value.transforms[i].translation.y;
        poses(2, i) = _value.transforms[i].translation.z;
        poses(3, i) = _value.transforms[i].rotation.w;
        poses(4, i) = _value.transforms[i].rotation.x;
        poses(5, i) = _value.transforms[i].rotation.y;
        poses(6, i) = _value.transforms[i].rotation.z;

        if(with_virtual)
        {
            poses(7, i) = _value.wrench[i].force.x;
            poses(8, i) = _value.wrench[i].force.y;
            poses(9, i) = _value.wrench[i].force.z;
        }
    }
    return poses;
}

DM KiteVisualizer::world2rviz(const DM &world_pose)
{
    DM pose = DM::zeros(world_pose.size());
    DM transform = DM::vertcat({0, 1, 0, 0});
    /** do only position transform if receive a point*/
    if(pose.size1() == 3)
    {
        DM q_pose = quat_multiply(transform, DM::vertcat({0, world_pose}) );
        DM tmp = quat_multiply(q_pose, quat_inverse(transform) );
        pose = tmp(Slice(1,4), 0);

        return pose;
    }

    /** position transform */
    DM q_pose = quat_multiply(transform, DM::vertcat({0, world_pose(Slice(0,3), 0) }) );
    DM tmp = quat_multiply(q_pose, quat_inverse(transform) );
    pose(Slice(0,3), 0) = tmp(Slice(1,4), 0);
    /** attitude transform */
    pose(Slice(3,7),0) = quat_multiply(quat_multiply(transform, world_pose(Slice(3,7), 0)), quat_inverse(transform));

    if(pose.size1() == 10)
    {
        q_pose = quat_multiply(transform, DM::vertcat({0, world_pose(Slice(7,10), 0) }) );
        tmp = quat_multiply(q_pose, quat_inverse(transform) );
        pose(Slice(7,10), 0) = tmp(Slice(1,4), 0);
    }

    return pose;
}

visualization_msgs::MarkerArray KiteVisualizer::getOptimalTrajectory()
{
    visualization_msgs::MarkerArray m_array;
    int array_size = optimal_trajectory.size2();
    if(array_size == 0)
        return m_array;

    DMVector split_traj = DM::horzsplit(optimal_trajectory, 1);
    visualization_msgs::Marker tmp_marker = m_kite_marker;
    int id_counter = m_kite_marker.id;

    DMVector::const_iterator it = split_traj.begin();
    while(std::distance<DMVector::const_iterator>(it, split_traj.end()) > 0)
    {
        DM state = (*it);
        //std::cout << state << "\n";
        id_counter += 1;

        std::vector<double> row = state.nonzeros();
        tmp_marker.pose.position.x = row[0];
        tmp_marker.pose.position.y = row[1];
        tmp_marker.pose.position.z = row[2];

        tmp_marker.pose.orientation.w = row[3];
        tmp_marker.pose.orientation.x = row[4];
        tmp_marker.pose.orientation.y = row[5];
        tmp_marker.pose.orientation.z = row[6];

        tmp_marker.id = id_counter;

        m_array.markers.push_back(tmp_marker);
        std::advance(it, 1);
    }
    return m_array;
}

visualization_msgs::Marker KiteVisualizer::getVirtualTrajectory()
{
    visualization_msgs::Marker m_array = virtual_state_markers;
    int array_size = optimal_trajectory.size2();
    if(array_size == 0)
        return m_array;

    DMVector split_traj = DM::horzsplit(optimal_trajectory, 1);
    m_array.points.clear();

    DMVector::const_iterator it = split_traj.begin();
    while( std::distance<DMVector::const_iterator>(it, split_traj.end()) > 0)
    {
        DM point = (*it)(Slice(7, 10));
        geometry_msgs::Point p;
        std::vector<double> coords = point.nonzeros();
        p.x = coords[0];
        p.y = coords[1];
        p.z = coords[2];

        m_array.points.push_back(p);
        std::advance(it, 1);
    }

    return m_array;
}

visualization_msgs::Marker KiteVisualizer::getTetherMarker()
{
    tether_marker.points[1] = m_kite_marker.pose.position;
    return tether_marker;
}

void KiteVisualizer::state2marker(const DM &kite_pose, visualization_msgs::Marker &_marker)
{
    std::vector<double> pose = kite_pose.nonzeros();
    _marker.pose.position.x = pose[0];
    _marker.pose.position.y = pose[1];
    _marker.pose.position.z = pose[2];

    _marker.pose.orientation.w = pose[3];
    _marker.pose.orientation.x = pose[4];
    _marker.pose.orientation.y = pose[5];
    _marker.pose.orientation.z = pose[6];
}


int main( int argc, char** argv )
{
    ros::init(argc, argv, "kite_visualization_node");
    ros::NodeHandle n;
    ros::Rate r(10);
    ros::Publisher marker_pub = n.advertise<visualization_msgs::Marker>("/visualization_marker", 10);
    ros::Publisher trajec_pub = n.advertise<visualization_msgs::MarkerArray>("/trajectory_marker", 10);

    KiteVisualizer kite(n);
    visualization_msgs::Marker gs_marker;
    visualization_msgs::Marker kite_marker;
    visualization_msgs::Marker path_marker;
    visualization_msgs::Marker tether_marker;

    MarkerProperties marker_props;
    marker_props.type  = visualization_msgs::Marker::CYLINDER;
    marker_props.scale = Scale(0.1, 0.1, 2.0);
    marker_props.color = Color(0.0, 1.0, 0.0, 1.0);
    marker_props.configureMarker(gs_marker);
    gs_marker.pose.position.z = 1.0;
    gs_marker.lifetime = ros::Duration();

    /** kite marker set up */
    marker_props.id = 1;
    marker_props.type = visualization_msgs::Marker::MESH_RESOURCE;
    marker_props.mesh_resource = "package://openkite/meshes/kite_rviz.dae";
    marker_props.scale = Scale(0.1, 0.1, 0.1);
    marker_props.configureMarker(kite_marker);
    kite_marker.lifetime = ros::Duration();

    /** path markers set up */
    marker_props.id = 3;
    marker_props.type = visualization_msgs::Marker::LINE_STRIP;
    marker_props.scale.x = 0.025;
    marker_props.color = Color(0, 1, 0, 0.5);
    marker_props.configureMarker(path_marker);
    path_marker.lifetime = ros::Duration();

    /** create a path */
    SX x = SX::sym("x");
    double h = M_PI / 6.0;
    double a = 0.2;
    double L = 5;
    SX theta = h + a * sin(2 * x);
    SX phi   = 4 * a * cos(x);
    SX Path  = SX::vertcat({theta, phi, L});

    Function path_fun = Function("path", {x}, {Path});
    for(double alpha = 0; alpha < 2 * M_PI; alpha += 0.2)
    {
        DMVector res = path_fun( DMVector{DM(alpha)} );
        DM sphere_point = res[0];
        DM point = kmath::spheric2cart<DM>(sphere_point(1), sphere_point(0), sphere_point(2));

        geometry_msgs::Point p;
        p.x = point.nonzeros()[0];
        p.y = point.nonzeros()[1];
        p.z = point.nonzeros()[2];
        path_marker.points.push_back(p);
    }

    bool published_once = false;

    while (ros::ok())
    {
        // Publish the marker
        while (marker_pub.getNumSubscribers() < 1)
        {
            if (!ros::ok())
            {
                return 0;
            }
            ROS_WARN_ONCE("Please create a subscriber to the marker");
            sleep(1);
        }

        if(kite.is_initialized())
        {
            kite_marker = kite.getKiteMarker();
            kite_marker.color.a = 1.0;
            marker_pub.publish(kite_marker);

            /** publish tether marker */
            tether_marker = kite.getTetherMarker();
            tether_marker.color.a = kite.tether_active ? 1.0 : 0.0;
            marker_pub.publish(tether_marker);
        }

        /** draw optimal trajectories */
        visualization_msgs::Marker virt_markers = kite.getVirtualTrajectory();
        visualization_msgs::MarkerArray opt_markers = kite.getOptimalTrajectory();

        //std::cout << "Kite : " << kite_marker << "\n";
        //std::cout << "Opt Traj : " << opt_markers << "\n";

        //if(!opt_markers.markers.empty())
        //    trajec_pub.publish(opt_markers);

        if(!virt_markers.points.empty())
            marker_pub.publish(virt_markers);

        if(!published_once)
        {
            //marker_pub.publish(gs_marker);
            marker_pub.publish(path_marker);
            published_once = true;
        }

        ros::spinOnce();
        r.sleep();
    }
}
