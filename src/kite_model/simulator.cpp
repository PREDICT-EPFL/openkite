#include "simulator.h"

using namespace casadi;

void Simulator::controlCallback(const openkite::aircraft_controls::ConstPtr &msg)
{
    // redo
    controls(0) = msg->thrust;
    controls(1) = msg->elevator;
    controls(2) = msg->ailerons;
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

    std::cout << "Control: " << controls << "\n";

    state = m_object->solve(state, controls, dt);
}

void Simulator::publish_state()
{
    std::vector<double> state_vec = state.nonzeros();

    msg_state.header.stamp = ros::Time::now();

    msg_state.twist[0].linear.x = state_vec[0];
    msg_state.twist[0].linear.y = state_vec[1];
    msg_state.twist[0].linear.z = 0.0;

    msg_state.twist[0].angular.x = 0.0; //state_vec[6];  // front wheel speed;
    msg_state.twist[0].angular.y = 0.0; //state_vec[7];  // rear wheel speed
    msg_state.twist[0].angular.z = state_vec[2];

    msg_state.transforms[0].translation.x = state_vec[3];
    msg_state.transforms[0].translation.y = state_vec[4];
    msg_state.transforms[0].translation.z = 0.0;

    DM angle = state_vec[5];
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
    KinematicRobot robot;

    int broadcast_state;
    n.param<int>("broadcast_state", broadcast_state, 1);

    /** create an integrator instance */
    double sim_rate;
    n.param<double>("simulation_rate", sim_rate, 100);
    /** cast to seconds and round to ms */
    double dt = (1/sim_rate);
    dt = std::roundf(dt * 1000) / 1000;

    Dict params({{"tf", dt}, {"tol", 1e-6}, {"method", CVODES}});
    Function ode = robot.getDynamics();
    ODESolver object(ode, params);

    Simulator simulator(object, n);
    ros::Rate loop_rate(100); /** 50 Hz */


    /** hacky control */
    /** read the telemetry log */
    DMVector telemetry;
    std::ifstream telemetry_file("/Users/plistov/EPFL/ROS/ros_catkin_ws/devel/lib/openkite/mpc_log2.txt", std::ios::in);

    /** load control data */
    if(!telemetry_file.fail())
    {
        while(!telemetry_file.eof())
        {
            std::vector<double> line;
            double entry;
            for (int i = 0; i < 26; ++i)
            {
                telemetry_file >> entry;
                line.push_back(entry);
            }
            telemetry.push_back(DM({line}).T());
        }
    }
    else
    {
        std::cout << "Could not open : telemetry file \n";
        telemetry_file.clear();
    }

    //DM ctl;// = telemetry[0](Slice(1, 4));
    //simulator.setControls(ctl);

    //int iter = 0;
    while (ros::ok())
    {
        /** read from the log and set controls */
        /**
        if(iter >= telemetry.size())
            iter = telemetry.size() - 1;
        ctl = telemetry[iter](Slice(1, 4));
        simulator.setControls(ctl);
        iter++;
        */
        //ctl = DM({0,0,0});
        //simulator.setControls(ctl);

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
