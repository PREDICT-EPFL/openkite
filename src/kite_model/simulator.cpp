#include "simulator.h"

using namespace casadi;

void Simulator::controlCallback(const openkite::aircraft_controls::ConstPtr &msg)
{
    controls[0] = msg->thrust;
    controls[1] = msg->elevator;
}

Simulator::Simulator(const ODESolver &object, const ros::NodeHandle &nh)
{
    m_object = std::make_shared<ODESolver>(object);
    m_nh     = std::make_shared<ros::NodeHandle>(nh);

    /** define dimensions first given solver object */
    controls = DM::zeros(m_object->dim_u());
    state    = DM::zeros(m_object->dim_x());

    /** initialize subscribers and publishers */
    int broadcast_state;
    m_nh->param<int>("broadcast_state", broadcast_state, 1);

    std::vector<double> initial_value;
    m_nh->getParam("init_state", initial_value);
    initialize(DM(initial_value));
    ROS_INFO_STREAM("Simulator initialized at: " << initial_value);

    //pose_pub  = m_nh->advertise<geometry_msgs::PoseStamped>("/kite_pose", 100);
    state_pub = m_nh->advertise<sensor_msgs::MultiDOFJointState>("/kite_state", 100);

    std::string control_topic = "/kite_controls";
    control_sub = m_nh->subscribe(control_topic, 100, &Simulator::controlCallback, this);

    msg_state.twist.resize(1);
    msg_state.transforms.resize(1);
    msg_state.joint_names.resize(1);
    msg_state.header.frame_id = "kite_sim";
    msg_state.joint_names[0] = "lox";
}

void Simulator::simulate()
{
    Dict p = m_object->getParams();
    double dt = p["tf"];
    //std::cout << "Received control signal: " << controls << "\n";
    //std::cout << "STATE_SIM: " << state << "\n";

    state = m_object->solve(state, controls, dt);
}

void Simulator::publish_state()
{
    std::vector<double> state_vec = state.nonzeros();

    msg_state.header.stamp = ros::Time::now();

    msg_state.twist[0].linear.x = 0.0;
    msg_state.twist[0].linear.y = 0.0;
    msg_state.twist[0].linear.z = 0.0;

    msg_state.twist[0].angular.x = 0.0;
    msg_state.twist[0].angular.y = 0.0;
    msg_state.twist[0].angular.z = 0.0;

    msg_state.transforms[0].translation.x = state_vec[0];
    msg_state.transforms[0].translation.y = state_vec[1];
    msg_state.transforms[0].translation.z = 0.0;

    DM angle = state_vec[2];
    DM q = DM::vertcat({cos(angle / 2), 0, 0, sin(angle / 2)});
    std::vector<double> q_vec = q.nonzeros();

    msg_state.transforms[0].rotation.w = q_vec[0];
    msg_state.transforms[0].rotation.x = 0.0;
    msg_state.transforms[0].rotation.y = 0.0;
    msg_state.transforms[0].rotation.z = q_vec[3];

    /** publish current state estimation */
    state_pub.publish(msg_state);
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "simulator");
    ros::NodeHandle n("~");

    /** create a kite object */
    MobileRobot robot;

    int broadcast_state;
    n.param<int>("broadcast_state", broadcast_state, 1);

    /** create an integrator instance */
    double sim_rate;
    n.param<double>("simulation_rate", sim_rate, 50);
    /** cast to seconds and round to ms */
    double dt = (1/sim_rate);
    dt = std::roundf(dt * 1000) / 1000;

    Dict params({{"tf", dt}, {"tol", 1e-6}, {"method", CVODES}});
    Function ode = robot.getDynamics();
    ODESolver object(ode, params);

    Simulator simulator(object, n);
    ros::Rate loop_rate(50); /** 50 Hz */

    while (ros::ok())
    {
        if(simulator.is_initialized())
        {
            simulator.simulate();
            if(broadcast_state)
                simulator.publish_state();
        }

        ros::spinOnce();
        loop_rate.sleep();
    }

    return 0;
}
