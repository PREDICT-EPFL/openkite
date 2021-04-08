#include "ekf_node.h"

using namespace casadi;

void KiteEKF_Node::filterCallback(const geometry_msgs::PoseStamped::ConstPtr &msg)
{
    /** update measurement */
    last_measurement_arrived = (*msg).header.stamp;
    //boost::unique_lock<boost::mutex> scoped_lock(m_mutex);
    /** @todo: revise sending mechanism instead */
    if ( !measurements.empty() )
    {
        double time_lapsed = (last_measurement_arrived - measurements.back().header.stamp).toSec();
        if( time_lapsed >= 0.01 )
        {
            measurements.push_back(*msg);
            //std::cout << "callback: measurement came at: " << last_measurement_arrived
            //          << " it was: " << convertToDM(*msg) << "\n";
        }
    }
    else
    {
        measurements.push_back(*msg);
        //std::cout << "callback: first measurement came at: " << last_measurement_arrived
        //          << " it was: " << convertToDM(*msg) << "\n";
    }
}

void KiteEKF_Node::controlCallback(const std_msgs::Int16MultiArray::ConstPtr &msg)
{
    double throttle = static_cast<double>(msg->data[0]);
    double elevator = static_cast<double>(msg->data[2]);
    double rudder = static_cast<double>(msg->data[3]);
    control = DM::vertcat({throttle, elevator, rudder});
}

KiteEKF_Node::KiteEKF_Node(const ros::NodeHandle &_nh, const Function &Integrator, const Function &Jacobian)
{
    /** @redo model initialization*/
    nh = std::make_shared<ros::NodeHandle>(_nh);
    measurements.set_capacity(2);

    kite_state.twist.resize(1);
    kite_state.transforms.resize(1);
    kite_state.joint_names.resize(1);
    kite_state.header.frame_id = "kite";
    kite_state.joint_names[0] = "lox";

    brf_offset = DMVector{-0.09, 0, -0.02};
    brf_rotation = DMVector{0.0, -1.0, 0.0, 0.0};

    /** create filter object */
    filter = std::make_shared<KiteEKF>(Integrator, Jacobian);

    /** initialize subscribers and publishers */
    state_pub = nh->advertise<sensor_msgs::MultiDOFJointState>("/kite_state", 100);

    /** @todo : parametrize topics */
    std::string pose_topic = "/optitrack_client/Kite/pose";
    std::string control_topic = "/chatter";

    pose_sub = nh->subscribe(pose_topic, 1000, &KiteEKF_Node::filterCallback, this);
    control_sub = nh->subscribe(control_topic, 1000, &KiteEKF_Node::controlCallback, this);

    m_initialized = false;
}

void KiteEKF_Node::initialize()
{
    /** initialize kite state based on 2 sequential pose measurements */
    if (measurements.size() == 2)
    {
        DM m_prev, m_new;
        double dt;

        /** @todo: work out proper mutex lock */
        //boost::unique_lock<boost::mutex> scoped_lock(m_mutex);
        m_prev = optitrack2world(convertToDM(measurements[0]));
        m_new = optitrack2world(convertToDM(measurements[1]));
        dt = measurements[1].header.stamp.toSec() - measurements[0].header.stamp.toSec();

        std::cout << "time: " << measurements[0].header.stamp << " meas_prev: " << m_prev << "\n";
        std::cout << "time: " << measurements[1].header.stamp << " meas_new: " << m_new << "\n";

        /** linear velocity estimation */
        DM rdot = ( m_new(Slice(0, 3), 0) - m_prev(Slice(0, 3), 0) ) / dt;
        SX att = m_prev( Slice(3, 7), 0);
        DM att_inv = kmath::quat_inverse(att);

        SX rdot_b = kmath::quat_multiply(att_inv, SX::vertcat({0, rdot}));
        DM v_body = kmath::quat_multiply(rdot_b, att);
        v_body = v_body(Slice(1, 4), 0);

        /** estimate angular rates */
        SX att2 = m_new( Slice(3, 7), 0);
        DM dq = kmath::quat_multiply(SX(att_inv), att2);
        DM identity = DM::vertcat({1,0,0,0});
        DM qw = (2.0 / dt) * (dq - identity);
        DM w_body = qw(Slice(1, 4), 0);

        std::cout << "Initialized at: " << v_body << " " << w_body << "\n";

        std::vector<double> v_brf = v_body.nonzeros();
        std::vector<double> w_brf = w_body.nonzeros();
        std::vector<double> pose_irf = m_new.nonzeros();

        kite_state.twist[0].linear.x = v_brf[0];
        kite_state.twist[0].linear.y = v_brf[1];
        kite_state.twist[0].linear.z = v_brf[2];

        kite_state.twist[0].angular.x = w_brf[3];
        kite_state.twist[0].angular.y = w_brf[4];
        kite_state.twist[0].angular.z = w_brf[5];

        kite_state.transforms[0].translation.x = pose_irf[0];
        kite_state.transforms[0].translation.y = pose_irf[1];
        kite_state.transforms[0].translation.z = pose_irf[2];

        kite_state.transforms[0].rotation.w = pose_irf[3];
        kite_state.transforms[0].rotation.x = pose_irf[4];
        kite_state.transforms[0].rotation.y = pose_irf[5];
        kite_state.transforms[0].rotation.z = pose_irf[6];

        /** update filter time */
        filter->setTime( measurements.back().header.stamp.toSec() );
        /** initialize filter */
        filter->setEstimation( DM::vertcat({v_body, w_body, m_new}) );

        m_initialized = true;
    }

}

DM KiteEKF_Node::convertToDM(const geometry_msgs::PoseStamped &_value)
{
    DM value = DM::zeros(7);
    value(0) = _value.pose.position.x;
    value(1) = _value.pose.position.y;
    value(2) = _value.pose.position.z;
    value(3) = _value.pose.orientation.w;
    value(4) = _value.pose.orientation.x;
    value(5) = _value.pose.orientation.y;
    value(6) = _value.pose.orientation.z;

    return value;
}

DM KiteEKF_Node::optitrack2world(const DM &opt_pose)
{
    DM position = opt_pose(Slice(0, 3), 0);
    DM attitude = opt_pose(Slice(3, 7), 0);

    DM tmp = kmath::quat_multiply(attitude, DM::vertcat({0, brf_offset}));
    DM offset_b = kmath::quat_multiply(tmp, kmath::quat_inverse(attitude));
    position += offset_b(Slice(1, 4), 0);

    /** transform to world ref frame */
    tmp = kmath::quat_multiply(brf_rotation, DM::vertcat({0, position}));
    DM world_pos = kmath::quat_multiply(tmp, kmath::quat_inverse(brf_rotation));
    world_pos = world_pos(Slice(1, 4), 0);

    tmp = kmath::quat_multiply(brf_rotation, attitude);
    DM world_att = kmath::quat_multiply(tmp, kmath::quat_inverse(brf_rotation));

    /** @todo : remove this hack */
    world_pos(2) = world_pos(2) + 2.77;

    return DM::vertcat({world_pos, world_att});
}

void KiteEKF_Node::estimate()
{
    //get observation
    DM observation = optitrack2world(convertToDM(measurements.back()));
    double tstamp = measurements.back().header.stamp.toSec();
    filter->estimate(observation, tstamp);

    /** update filter time */
    filter->setTime(tstamp);
}

void KiteEKF_Node::publish()
{
    /** pack estimation to ROS message */
    DM correction = DM({0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        -0.09, -0.1247, -0.0418, 0.0, 0.0, 0.0, 0.0}) +
            DM({0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    -0.02, 0.01, -0.0418, 0.0, 0.0, 0.0, 0.0});
    DM corrected_estim = filter->getEstimation() - correction;
    //std::cout << "Estimated distance: " << DM::norm_2(corrected_estim(Slice(6,9))) << "\n";
    std::vector<double> state_vec = corrected_estim.nonzeros();

    kite_state.header.stamp = ros::Time::now();

    kite_state.twist[0].linear.x = state_vec[0];
    kite_state.twist[0].linear.y = state_vec[1];
    kite_state.twist[0].linear.z = state_vec[2];

    kite_state.twist[0].angular.x = state_vec[3];
    kite_state.twist[0].angular.y = state_vec[4];
    kite_state.twist[0].angular.z = state_vec[5];

    kite_state.transforms[0].translation.x = state_vec[6];
    kite_state.transforms[0].translation.y = state_vec[7];
    kite_state.transforms[0].translation.z = state_vec[8];

    kite_state.transforms[0].rotation.w = state_vec[9];
    kite_state.transforms[0].rotation.x = state_vec[10];
    kite_state.transforms[0].rotation.y = state_vec[11];
    kite_state.transforms[0].rotation.z = state_vec[12];

    /** publish current state estimation */
    state_pub.publish(kite_state);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "ekf_node");
    ros::NodeHandle n("~");

    unsigned tries = 0;
    unsigned locks = 0;

    /** create a kite object */
    std::string kite_params_file;
    n.param<std::string>("kite_params", kite_params_file, "./umx_radian.yaml");
    //std::cout << kite_params_file << "\n";
    //KiteProperties kite_props = kite_utils::LoadProperties(kite_params_file);
    AlgorithmProperties algo_props;
    algo_props.Integrator = RK4;
    algo_props.sampling_time = 0.02;
    //KiteDynamics kite = KiteDynamics(kite_props, algo_props);

    /** create a rigid body topic */
    RigidBodyKinematics rigid_body = RigidBodyKinematics(algo_props);

    Function fIntegrator = rigid_body.getNumericIntegrator();
    Function fJacobian = rigid_body.getNumericJacobian();

    /** create an EKF instance */
    KiteEKF_Node filter(n, fIntegrator, fJacobian);
    ros::Rate loop_rate(50); /** 50 Hz */

    while (ros::ok())
    {
        //std::cout << "TRIES: " << tries << " LOCKS: " << locks << "\n";

        /** Some complicated logic here please */
        if (!filter.initialized())
        {   std::cout << "CANNOT INITIALIZE \n";
            //boost::unique_lock<boost::mutex> scoped_lock(filter.m_mutex);
            filter.initialize();
        }
        else
        {
            /** perform estimation */
            //tries++;
            //boost::unique_lock<boost::mutex> scoped_lock(filter.m_mutex);
            ///locks++;
            filter.estimate();
            filter.publish();
        }

        ros::spinOnce();
        loop_rate.sleep();
    }

    return 0;
}
