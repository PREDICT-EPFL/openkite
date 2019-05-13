#include "nmpf_node.hpp"
#include "ekf_node.h"


using namespace casadi;

void KiteNMPF_Node::filterCallback(const sensor_msgs::MultiDOFJointState::ConstPtr &msg)
{
    //boost::unique_lock<boost::mutex> scoped_lock(m_mutex, boost::try_to_lock);
    DM estimation = convertToDM(*msg);
    /** prevent outliers*/
    std::vector<double> estim = estimation.nonzeros();
    if((std::fabs(estim[3]) >= 4) || (std::fabs(estim[4]) >= 7) || (std::fabs(estim[5]) >= 4))
    {
    }
    else
    {
        kite_state = estimation;
    }
    if(!is_initialized())
        initialize();
}

KiteNMPF_Node::KiteNMPF_Node(const ros::NodeHandle &_nh, const KiteProperties &kite_props,
                                                         const AlgorithmProperties &algo_props )
{
    std::shared_ptr<KiteDynamics> kite = std::make_shared<KiteDynamics>(kite_props, algo_props);

    /** @badcode : parametrize path outside NMPF node  constructor*/
    SX x = SX::sym("x");
    double radius   = 2.65;
    double altitude = 0.00;
    SX Path = SX::vertcat(SXVector{radius * cos(x), radius * sin(x), altitude});
    /** rotate path */
    SX q_rot = SX::vertcat({cos(M_PI / 8), 0, sin(M_PI / 8), 0});
    SX q_rot_inv = kmath::quat_inverse(q_rot);
    SX qP_tmp = kmath::quat_multiply(q_rot_inv, SX::vertcat({0, Path}));
    SX qP_q = kmath::quat_multiply(qP_tmp, q_rot);
    Path = qP_q(Slice(1,4), 0);
    Function path = Function("path", {x}, {Path});

    controller = std::make_shared<KiteNMPF>(kite, path);
    nh = std::make_shared<ros::NodeHandle>(_nh);
    /** set control constraints */
    double angle_sat = kmath::deg2rad(7.0);
    DM lbu = DM::vertcat({0.1, -angle_sat, -angle_sat, -5});
    DM ubu = DM::vertcat({0.15, angle_sat, angle_sat, 5});

    /** scaling matrices */
    DM ScaleX  = DM::diag(DM({0.1, 1/3.0, 1/3.0, 1/2.0, 1/5.0, 1/2.0, 1/3.0, 1/3.0, 1/3.0, 1.0, 1.0, 1.0, 1.0, 1/6.28, 1/6.28}));
    DM ScaleU = DM::diag(DM({1/0.15, 1/0.2618, 1/0.2618, 1/5.0}));
    controller->setControlScaling(ScaleU);
    controller->setStateScaling(ScaleX);

    controller->setLBU(lbu);
    controller->setUBU(ubu);

    /** set variable constraints */
    DM lbx = DM::vertcat({2.0, -DM::inf(1), -DM::inf(1), -4 * M_PI, -4 * M_PI, -4 * M_PI, -DM::inf(1), -DM::inf(1), -DM::inf(1),
                         -1.01, -1.01, -1.01, -1.01, -DM::inf(1), -DM::inf(1)});

    DM ubx = DM::vertcat({DM::inf(1), DM::inf(1), DM::inf(1), 4 * M_PI, 4 * M_PI, 4 * M_PI, DM::inf(1), DM::inf(1), DM::inf(1),
                          1.01, 1.01, 1.01, 1.01, DM::inf(1), DM::inf(1)});

    controller->setLBX(lbx);
    controller->setUBX(ubx);

    DM vel_ref = 4.0;
    controller->setReferenceVelocity(vel_ref);

    /** create NLP */
    controller->createNLP();

    /** create solver for delay compensation */
    nh->param<double>("delay", transport_delay, 0.1);

    Dict opts;
    opts["tf"]         = transport_delay;
    opts["poly_order"] = 2;
    opts["tol"]        = 1e-4;
    opts["method"]     = IntType::CVODES;
    Function system    = kite->getNumericDynamics();
    //Function system = controller->getAugDynamics();
    solver = std::make_shared<ODESolver>(system, opts);

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
    control    = opt_ctl(Slice(0,3), opt_ctl.size2() - 1);
    std::vector<double> controls = control.get_nonzeros();
    openkite::aircraft_controls control_msg;

    if(!controls.empty())
    {
        control_msg.header.stamp = ros::Time::now();
        control_msg.thrust = controls[0];
        control_msg.elevator = controls[1];
        control_msg.rudder = controls[2];

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
        opt_msg.twist[idx].linear.x = row[0];
        opt_msg.twist[idx].linear.y = row[1];
        opt_msg.twist[idx].linear.z = row[2];

        opt_msg.twist[idx].angular.x = row[3];
        opt_msg.twist[idx].angular.y = row[4];
        opt_msg.twist[idx].angular.z = row[5];

        opt_msg.transforms[idx].translation.x = row[6];
        opt_msg.transforms[idx].translation.y = row[7];
        opt_msg.transforms[idx].translation.z = row[8];

        opt_msg.transforms[idx].rotation.w = row[9];
        opt_msg.transforms[idx].rotation.x = row[10];
        opt_msg.transforms[idx].rotation.y = row[11];
        opt_msg.transforms[idx].rotation.z = row[12];

        /** virtual state */
        Function path = controller->getPathFunction();
        DM point  = path(DMVector{row[13]})[0];
        std::vector<double> virt_point = point.nonzeros();

        opt_msg.wrench[idx].force.x = virt_point[0];
        opt_msg.wrench[idx].force.y = virt_point[1];
        opt_msg.wrench[idx].force.z = virt_point[2];
        idx++;
    }

    traj_pub.publish(opt_msg);
}

/** publish diagnostic info */
void KiteNMPF_Node::publish_mpc_diagnostic()
{
    openkite::mpc_diagnostic diag_msg;
    diag_msg.header.stamp = ros::Time::now();

    diag_msg.pos_error    = controller->getPathError();
    diag_msg.comp_time_ms = comp_time_ms * 1000;
    diag_msg.virt_state   = controller->getVirtState();
    diag_msg.vel_error    = controller->getVelocityError();

    /** dummy output */
    diag_msg.cost         = 0;
    diagnostic_pub.publish(diag_msg);
}

void KiteNMPF_Node::compute_control()
{
    /** augment state with pseudo state*/
    DM augmented_state;
    DM opt_traj = controller->getOptimalTrajetory();

    /** make local copy */
    DM local_copy = kite_state;

    if(!opt_traj.is_empty())
    {
        /** transport delay compensation */
        DM predicted_state = solver->solve(local_copy, control, transport_delay);
        //std::cout << "virtual state : " << opt_traj << "\n";
        augmented_state = DM::vertcat({predicted_state, opt_traj(Slice(13, opt_traj.size1()), opt_traj.size2() - 3)});

        Dict stats = controller->getStats();
        std::string solve_status = static_cast<std::string>(stats["return_status"]);

//        Dict stats = controller->getStats();
//        std::string solve_status = static_cast<std::string>(stats["return_status"]);
//        if(solve_status.compare("Solve_Succeeded") != 0)
//        {
//            std::cout << "REINIT! \n";
//            augmented_state(13) = controller->findClosestPointOnPath(kite_state(Slice(6,9)), augmented_state(13));
//        }
    }
    else
    {
        DM closest_point = controller->findClosestPointOnPath(kite_state(Slice(6,9)));
        std::cout << "Initialiaztion: " << closest_point << "\n";
        augmented_state = DM::vertcat({kite_state, closest_point, 0});
    }

    /** @badcode : dirty hack : zero-speed workaround */
    DM minimal_speed = DM(2.1);
    if (augmented_state(0, 0).nonzeros()[0] < minimal_speed.nonzeros()[0])
        augmented_state(0, 0) = minimal_speed;
    /** compute control */
    controller->computeControl(augmented_state);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "nmpf_node");
    ros::NodeHandle n;

    /** create a kite object */
    std::string kite_params_file;
    n.param<std::string>("kite_params", kite_params_file, "umx_radian2.yaml");
    std::cout << kite_params_file << "\n";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_params_file);
    AlgorithmProperties algo_props;
    algo_props.Integrator = RK4;
    algo_props.sampling_time = 0.02;

    int broadcast_trajectory;
    n.param<int>("broadcast_trajectory", broadcast_trajectory, 1);

    /** create a NMPF instance */
    KiteNMPF_Node tracker(n, kite_props, algo_props);
    ros::Rate loop_rate(14); /** 18 Hz */

    while (ros::ok())
    {
        ros::spinOnce();

        if(tracker.is_initialized())
        {
            double start = ros::Time::now().toSec();
            tracker.compute_control();
            double finish = ros::Time::now().toSec();
            tracker.publish();
            tracker.publish_mpc_diagnostic();
            tracker.comp_time_ms = finish - start;
            std::cout << "Control computational delay: " << finish - start << "\n";

            if(broadcast_trajectory)
                tracker.publish_trajectory();
        }
        else
        {
            loop_rate.sleep();
        }

    }

    return 0;
}
