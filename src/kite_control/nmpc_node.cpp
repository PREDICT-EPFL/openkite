#include "nmpc_node.hpp"
#include "ekf_node.h"


void KiteNMPF_Node::filterCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg)
{

    DM estimation = convertToDM(*msg);
    kite_state = casadi::DM::zeros(6);
    kite_state(0) = estimation(0); // vx
    kite_state(1) = estimation(1); // vy
    kite_state(2) = estimation(5); // omega

    kite_state(3) = estimation(6); // x
    kite_state(4) = estimation(7); // y
    kite_state(5) = 2 * asin(estimation(12)); // theta

    //kite_state[6] = estimation[3];  // omega_fw
    //kite_state[7] = estimation[4];  // omega_rw

    if(!is_initialized())
        initialize();

    // kinematic robot
    /**
    DM estimation = convertToDM(*msg);
    kite_state = casadi::DM::zeros(3);

    kite_state[0] = estimation[6]; // x
    kite_state[1] = estimation[7]; // y
    kite_state[2] = 2 * asin(estimation[12]); // theta

    if(!is_initialized())
        initialize();
        */
}


KiteNMPF_Node::KiteNMPF_Node(const ros::NodeHandle &_nh)
{
    /** create nmpf controller here */
    //const int dimx = 6;
    //const int dimu = 3;
    double tf = 2.0;

    DMDict options;
    options["mpc.scaling"] = 1;
    options["mpc.scale_x"] = DM::diag(DM({1,1,1,1,1,1}));
    options["mpc.scale_u"] = DM::diag(DM({1.0,1.0/4000,1.0/4000}));

    controller = std::make_shared<NMPC>(casadi::DM::zeros(3),tf, options);

    /** set state and control constraints */
    DM lbu = DM::vertcat({-0.8, 0.0, -2000});
    DM ubu = DM::vertcat({ 0.8, 0.0,  2000});
    controller->setLBU(lbu);
    controller->setUBU(ubu);

    /** set state/virtual state constraints */
    DM lbx = DM::vertcat({-1.5, -DM::inf(), -DM::inf(), -DM::inf(), -DM::inf(), -DM::inf()});
    DM ubx = DM::vertcat({1.5,  DM::inf(),  DM::inf(),  DM::inf(),  DM::inf(),  DM::inf()});
    controller->setLBX(lbx);
    controller->setUBX(ubx);

    nh = std::make_shared<ros::NodeHandle>(_nh);

    /** initialize subscribers and publishers */
    control_pub    = nh->advertise<openkite::aircraft_controls>("kite_controls", 100);
    traj_pub       = nh->advertise<sensor_msgs::MultiDOFJointState>("/opt_traj", 10);
    diagnostic_pub = nh->advertise<openkite::mpc_diagnostic>("/mpc_diagnostic", 10);

    std::string state_topic = "/kite_state";
    state_sub = nh->subscribe(state_topic, 100, &KiteNMPF_Node::filterCallback, this);

    m_initialized = false;
    comp_time_ms  = 0.0;
}


DM KiteNMPF_Node::convertToDM(const sensor_msgs::MultiDOFJointState &_value)
{
    DM value = DM::zeros(13);
    value(0) = _value.twist.back().linear.x;
    value(1) = _value.twist.back().linear.y;
    value(2) = _value.twist.back().linear.z;
    value(3) = _value.twist.back().angular.x;
    value(4) = _value.twist.back().angular.y;
    value(5) = _value.twist.back().angular.z;
    value(6) = _value.transforms.back().translation.x;
    value(7) = _value.transforms.back().translation.y;
    value(8) = _value.transforms.back().translation.z;
    value(9) = _value.transforms.back().rotation.w;
    value(10) = _value.transforms.back().rotation.x;
    value(11) = _value.transforms.back().rotation.y;
    value(12) = _value.transforms.back().rotation.z;

    return value;
}


void KiteNMPF_Node::publish()
{
    /** pack estimation to ROS message */
    DM opt_ctl = controller->getOptimalControl();
    control    = opt_ctl(Slice(0, 3), opt_ctl.size2() - 2);
    std::vector<double> controls = control.get_nonzeros();
    openkite::aircraft_controls control_msg;

    if(!controls.empty())
    {
        control_msg.header.stamp = ros::Time::now();
        control_msg.thrust = controls[0];
        control_msg.elevator = controls[1];
        control_msg.ailerons = controls[2];

        control_msg.flaps = 0.0;
        control_msg.rudder = 0.0;

        /** publish current control */
        control_pub.publish(control_msg);
    }
}

void KiteNMPF_Node::publish_mpc_diagnostic()
{
    openkite::mpc_diagnostic diag_msg;
    diag_msg.header.stamp = ros::Time::now();

    diag_msg.pos_error    = 0.0;
    diag_msg.comp_time_ms = comp_time_ms * 1000;
    diag_msg.virt_state   = 0.0;
    diag_msg.vel_error    = 0.0;

    /** dummy output */
    diag_msg.cost         = 0;
    diagnostic_pub.publish(diag_msg);
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

        opt_msg.transforms[idx].translation.x = row[3];
        opt_msg.transforms[idx].translation.y = row[4];
        opt_msg.transforms[idx].translation.z = 0.0;

        /*** attitude */
        DM angle = row[5];
        DM q = DM::vertcat({cos(angle / 2), 0, 0, sin(angle / 2)});
        std::vector<double> q_vec = q.nonzeros();

        opt_msg.transforms[idx].rotation.w = q_vec[0];
        opt_msg.transforms[idx].rotation.x = 0.0;
        opt_msg.transforms[idx].rotation.y = 0.0;
        opt_msg.transforms[idx].rotation.z = q_vec[3];

        opt_msg.wrench[idx].force.x = 0.0;
        opt_msg.wrench[idx].force.y = 0.0;
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

    if(!opt_traj.is_empty())
    {
        augmented_state = kite_state;

        Dict stats = controller->getStats();
        std::string solve_status = static_cast<std::string>(stats["return_status"]);
    }
    else
    {
        augmented_state = kite_state;
        //augmented_state = DM({0,0,0,-5,0,0});
    }
    /** compute control */
    controller->computeControl(augmented_state);
    //DM opt_trajectory = controller->getOptimalTrajetory();
    //DM opt_control    = controller->getOptimalControl();

    //std::cout << "Optimal Trajectory: \n" << opt_trajectory << "\n";
    //std::cout << "Optimal Control: \n" << opt_control << "\n";
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "nmpc_node");
    ros::NodeHandle n;

    int broadcast_trajectory;
    n.param<int>("broadcast_trajectory", broadcast_trajectory, 1);

    /** create a NMPF instance */
    KiteNMPF_Node tracker(n);
    ros::Rate loop_rate(50); /** 50 Hz */


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
            tracker.publish_mpc_diagnostic();
            std::cout << "Control computational delay: " << finish - start << "\n";

            if(broadcast_trajectory)
                tracker.publish_trajectory();
        }

        loop_rate.sleep();
    }

    return 0;
}
