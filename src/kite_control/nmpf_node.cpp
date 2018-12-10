#include "nmpf_node.hpp"
#include "ekf_node.h"


void KiteNMPF_Node::filterCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg)
{
    DM estimation = convertToDM(*msg);
    kite_state = casadi::DM::ones(3);
    kite_state[0] = estimation[6]; // x
    kite_state[1] = estimation[7]; // y
    kite_state[2] = 2 * asin(estimation[12]); // theta

    if(!is_initialized())
        initialize();
}


KiteNMPF_Node::KiteNMPF_Node(const ros::NodeHandle &_nh)
{
    /** create nmpf controller here */
    const int dimx = 3;
    const int dimu = 2;
    double tf = 2.0;
    controller = std::make_shared<NMPF>(tf);

    /** set state and control constraints */
    DM lbu = DM::vertcat({0.0, -0.7, -1.0});
    DM ubu = DM::vertcat({1.0, 0.7, 1.0});
    controller->setLBU(lbu);
    controller->setUBU(ubu);

    /** set state/virtual state constraints */
    DM lbx = DM::vertcat({-DM::inf(), -DM::inf(), -DM::inf(), 0.0, 0.0});
    DM ubx = DM::vertcat({DM::inf(), DM::inf(), DM::inf(), 11.3, DM::inf()});
    controller->setLBX(lbx);
    controller->setUBX(ubx);

    nh = std::make_shared<ros::NodeHandle>(_nh);

    DM vel_ref = 1.0;
    controller->setReferenceVelocity(vel_ref);

    /** initialize subscribers and publishers */
    control_pub    = nh->advertise<openkite::aircraft_controls>("kite_controls", 100);
    traj_pub       = nh->advertise<sensor_msgs::MultiDOFJointState>("/opt_traj", 10);
    diagnostic_pub = nh->advertise<openkite::mpc_diagnostic>("/mpc_diagnostic", 10);

    std::string state_topic = "/kite_state";
    state_sub = nh->subscribe(state_topic, 100, &KiteNMPF_Node::filterCallback, this);

    m_initialized = false;
    comp_time_ms = 0.0;
}


DM KiteNMPF_Node::convertToDM(const sensor_msgs::MultiDOFJointState &_value)
{
    DM value = DM::zeros(13);
    value[0] = _value.twist.back().linear.x;
    value[1] = _value.twist.back().linear.y;
    value[2] = _value.twist.back().linear.z;
    value[3] = _value.twist.back().angular.x;
    value[4] = _value.twist.back().angular.y;
    value[5] = _value.twist.back().angular.z;
    value[6] = _value.transforms.back().translation.x;
    value[7] = _value.transforms.back().translation.y;
    value[8] = _value.transforms.back().translation.z;
    value[9] = _value.transforms.back().rotation.w;
    value[10] = _value.transforms.back().rotation.x;
    value[11] = _value.transforms.back().rotation.y;
    value[12] = _value.transforms.back().rotation.z;

    return value;
}


void KiteNMPF_Node::publish()
{
    /** pack estimation to ROS message */
    DM opt_ctl = controller->getOptimalControl();
    control    = opt_ctl(Slice(0,2), opt_ctl.size2() - 1);
    std::vector<double> controls = control.get_nonzeros();
    openkite::aircraft_controls control_msg;

    if(!controls.empty())
    {
        control_msg.header.stamp = ros::Time::now();
        control_msg.thrust = controls[0];
        control_msg.elevator = controls[1];
        control_msg.ailerons = 0.0;
        control_msg.flaps = 0.0;
        control_msg.rudder = 0.0;

        /** publish current control */
        control_pub.publish(control_msg);
    }
}

void KiteNMPF_Node::publish_trajectory()
{
    DM opt_trajectory = controller->getOptimalTrajetory();

    int array_size = opt_trajectory.size2();
    if(array_size == 0)
        return;

    DMVector split_traj = DM::horzsplit(opt_trajectory, 1);
    sensor_msgs::MultiDOFJointState opt_msg;
    opt_msg.twist.resize(array_size);
    opt_msg.transforms.resize(array_size);
    opt_msg.joint_names.resize(array_size);
    opt_msg.wrench.resize(array_size);
    opt_msg.header.frame_id = "optimal_trajectory";
    int idx = 0;

    for(DMVector::const_iterator it = split_traj.begin(); it != split_traj.end(); std::advance(it, 1))
    {
        std::vector<double> row = (*it).nonzeros();
        opt_msg.twist[idx].linear.x = 0.0;
        opt_msg.twist[idx].linear.y = 0.0;
        opt_msg.twist[idx].linear.z = 0.0;

        opt_msg.twist[idx].angular.x = 0.0;
        opt_msg.twist[idx].angular.y = 0.0;
        opt_msg.twist[idx].angular.z = 0.0;

        opt_msg.transforms[idx].translation.x = row[0];
        opt_msg.transforms[idx].translation.y = row[1];
        opt_msg.transforms[idx].translation.z = 0.0;

        /*** attitude */
        DM angle = row[2];
        DM q = DM::vertcat({cos(angle / 2), 0, 0, sin(angle / 2)});
        std::vector<double> q_vec = q.nonzeros();

        opt_msg.transforms[idx].rotation.w = q_vec[0];
        opt_msg.transforms[idx].rotation.x = 0.0;
        opt_msg.transforms[idx].rotation.y = 0.0;
        opt_msg.transforms[idx].rotation.z = q_vec[3];

        /** virtual state */
        Function path = controller->getPathFunction();
        DM point  = path(DMVector{row[3]})[0];
        std::vector<double> virt_point = point.nonzeros();

        opt_msg.wrench[idx].force.x = virt_point[0];
        opt_msg.wrench[idx].force.y = virt_point[1];
        opt_msg.wrench[idx].force.z = 0.0;
        idx++;
    }

    traj_pub.publish(opt_msg);
}


void KiteNMPF_Node::compute_control()
{
    /** augment state with pseudo state*/
    DM augmented_state;
    DM opt_traj = controller->getOptimalTrajetory();

    //std::cout << "STATE: " << kite_state << "\n";

    if(!opt_traj.is_empty())
    {
        augmented_state = DM::vertcat({kite_state, opt_traj(Slice(3, opt_traj.size1()), opt_traj.size2() - 3)});

        Dict stats = controller->getStats();
        std::string solve_status = static_cast<std::string>(stats["return_status"]);
    }
    else
    {
        DM closest_point = controller->findClosestPointOnPath(kite_state(Slice(0,2)));
        if(closest_point.nonzeros()[0] <= 0)
            closest_point = 0.0001;
        std::cout << "Initialization: " << closest_point << "\n";
        augmented_state = DM::vertcat({kite_state, closest_point, 0});
    }
    /** compute control */
    controller->computeControl(augmented_state);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "nmpf_node");
    ros::NodeHandle n;

    int broadcast_trajectory;
    n.param<int>("broadcast_trajectory", broadcast_trajectory, 1);

    // sleep
    //ros::Duration(10.0).sleep();

    /** create a NMPF instance */
    KiteNMPF_Node tracker(n);
    ros::Rate loop_rate(20); /** 18 Hz */

    while (ros::ok())
    {
        ros::spinOnce();

        if(tracker.is_initialized())
        {
            double start = ros::Time::now().toSec();
            tracker.compute_control();
            double finish = ros::Time::now().toSec();
            tracker.publish();
            tracker.comp_time_ms = finish - start;
            std::cout << "Control computational delay: " << finish - start << "\n";

            if(broadcast_trajectory)
                tracker.publish_trajectory();
        }

        loop_rate.sleep();
    }

    return 0;
}
