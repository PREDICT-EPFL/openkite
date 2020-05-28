#include "simulator.h"

using namespace casadi;

//void Simulator::controlCallback(const openkite::aircraft_controls::ConstPtr &msg)
//{
//    controls(0) = msg->thrust;
//    controls(1) = msg->elevator;
//    controls(2) = msg->rudder;
//    controls(3) = msg->ailerons;
//}

void Simulator::controlCallback(const sensor_msgs::JoyConstPtr &msg) {
    controls(0) = msg->axes[0]; // Thrust
    controls(1) = msg->axes[1]; // Elevator
    controls(2) = msg->axes[2]; // Rudder
    controls(3) = msg->axes[3]; // Aileron
}

Simulator::Simulator(const ODESolver &odeSolver, const ros::NodeHandle &nh) {
    m_odeSolver = std::make_shared<ODESolver>(odeSolver);
    m_nh = std::make_shared<ros::NodeHandle>(nh);

    /** define dimensions first given solver object */
    controls = DM::zeros(m_odeSolver->dim_u());
    state = DM::zeros(m_odeSolver->dim_x());

    /** initialize subscribers and publishers */
    int broadcast_state;
    m_nh->param<int>("broadcast_state", broadcast_state, 1);

    /** initialize subscribers and publishers */
    m_nh->param<bool>("simulate_tether", sim_tether, false);

    std::vector<double> initial_value;
    m_nh->getParam("init_state", initial_value);
    initial_value.emplace_back(0);
    initial_value.emplace_back(0);
    initial_value.emplace_back(0);
    initial_value.emplace_back(0);
    initialize(DM(initial_value));
    ROS_INFO_STREAM("Simulator initialized at: " << initial_value);

    //pose_pub  = m_nh->advertise<geometry_msgs::PoseStamped>("/sim/kite_pose", 100);
    state_pub = m_nh->advertise<sensor_msgs::MultiDOFJointState>("/sim/kite_state", 1);
    tether_pub = m_nh->advertise<geometry_msgs::Vector3Stamped>("/sim/tether_force", 1);

    control_sub = m_nh->subscribe("/sim/set/kite_controls", 100, &Simulator::controlCallback, this);

    msg_state.twist.resize(2);
    msg_state.transforms.resize(1);
    msg_state.wrench.resize(2);

    msg_state.header.frame_id = "kite_sim";
}

void Simulator::simulate() {
    Dict p = m_odeSolver->getParams();
    double dt = p["tf"];

    /* Make sure thrust is not negative */
    if (static_cast<double>(controls(0)) < 0.0)
        controls(0) = 0.0;

    /* Get pitot airspeed */
    DM airspeed_evaluated = m_NumericAirspeedMeas(DMVector{state})[0];
    Va_pitot = airspeed_evaluated.nonzeros()[0];

    /* Get aero values (Va, alpha, beta) */
    DM aeroValues_evaluated = m_NumericAeroValues(DMVector{state})[0];
    std::vector<double> aeroValues_evaluated_vect = aeroValues_evaluated.nonzeros();
    Va = aeroValues_evaluated_vect[0];
    alpha = aeroValues_evaluated_vect[1];
    beta = aeroValues_evaluated_vect[2];

    /* Get specific nongravitational force before solving (altering) the state */
    DM specNongravForce_evaluated = m_NumericSpecNongravForce(DMVector{state, controls})[0];
    std::vector<double> specNongravForce_evaluated_vect = specNongravForce_evaluated.nonzeros();
    specNongravForce = specNongravForce_evaluated_vect;

    if (sim_tether) {
        /* Get tether force vector */
        DM specTethForce_evaluated = m_NumericSpecTethForce(DMVector{state, controls})[0];
        std::vector<double> specTethForce_evaluated_vect = specTethForce_evaluated.nonzeros();
        specTethForce = specTethForce_evaluated_vect;
    }

    state = m_odeSolver->solve(state, controls, dt);
    //std::cout << "length : " << DM::norm_2(state(Slice(6,9))) << "\n";
    //std::cout << "State: " << state << "\n";
}

void Simulator::publish_state(const ros::Time &sim_time) {
    /** State message */
    std::vector<double> state_vec = state.nonzeros();

    msg_state.header.stamp = sim_time;

    msg_state.twist[0].linear.x = state_vec[0];
    msg_state.twist[0].linear.y = state_vec[1];
    msg_state.twist[0].linear.z = state_vec[2];

    msg_state.twist[0].angular.x = state_vec[3];
    msg_state.twist[0].angular.y = state_vec[4];
    msg_state.twist[0].angular.z = state_vec[5];

    msg_state.transforms[0].translation.x = state_vec[6];
    msg_state.transforms[0].translation.y = state_vec[7];
    msg_state.transforms[0].translation.z = state_vec[8];

    msg_state.transforms[0].rotation.w = state_vec[9];
    msg_state.transforms[0].rotation.x = state_vec[10];
    msg_state.transforms[0].rotation.y = state_vec[11];
    msg_state.transforms[0].rotation.z = state_vec[12];

    /* Using wrench 0 as for acceleration, as measured by IMU (spec nongravitational forces) */
    msg_state.wrench[0].force.x = specNongravForce[0];
    msg_state.wrench[0].force.y = specNongravForce[1];
    msg_state.wrench[0].force.z = specNongravForce[2];

    /* Using wrench 1 as for acceleration */
    msg_state.wrench[1].force.x = specTethForce[0];
    msg_state.wrench[1].force.y = specTethForce[1];
    msg_state.wrench[1].force.z = specTethForce[2];

    /* Aero values */
    msg_state.twist[1].linear.x = Va;
    msg_state.twist[1].linear.y = beta;
    msg_state.twist[1].linear.z = alpha;
    msg_state.twist[1].angular.x = Va_pitot;

    state_pub.publish(msg_state);
}

void Simulator::publish_pose(const ros::Time &sim_time) {
    DM pose = getPose();
    std::vector<double> pose_vec = pose.nonzeros();

    geometry_msgs::PoseStamped msg_pose;
    msg_pose.header.stamp = sim_time;

    msg_pose.pose.position.x = pose_vec[0];
    msg_pose.pose.position.y = pose_vec[1];
    msg_pose.pose.position.z = pose_vec[2];

    msg_pose.pose.orientation.w = pose_vec[3];
    msg_pose.pose.orientation.x = pose_vec[4];
    msg_pose.pose.orientation.y = pose_vec[5];
    msg_pose.pose.orientation.z = pose_vec[6];

    pose_pub.publish(msg_pose);
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "simulator");
    ros::NodeHandle n("~");

    /** create a kite object */
    std::string kite_params_file;
    n.param<std::string>("kiteparams_path", kite_params_file, "_kiteparams_.yaml");
    std::cout << "Simulator: Using kite parameters from: " << kite_params_file << "\n";
    KiteProperties kite_props = kite_utils::LoadProperties(kite_params_file);

    double windFrom_deg{180};
    n.param<double>("windFrom_deg", windFrom_deg, 180.0);
    kite_props.atmosphere.WindFrom = windFrom_deg * M_PI / 180.0;
    n.param<double>("windSpeed", kite_props.atmosphere.WindSpeed, 0.0);
    kite_props.atmosphere.airDensity = 1.1589; // Standard atmosphere at 468 meters

    std::cout << "Simulator: Wind from " << windFrom_deg << " deg at " << kite_props.atmosphere.WindSpeed << " m/s.\n";

    bool simulate_tether;
    n.param<bool>("simulate_tether", simulate_tether, false);
    if (simulate_tether)
        std::cout << "[Simulator]: Tether ON\n";

    AlgorithmProperties algo_props;
    algo_props.Integrator = RK4;
    algo_props.sampling_time = 0.02;
    KiteDynamics kite = KiteDynamics(kite_props, algo_props, simulate_tether);
//    MinimalKiteDynamics kite = MinimalKiteDynamics(kite_props, algo_props);

    int broadcast_state;
    n.param<int>("broadcast_state", broadcast_state, 1);

    /** create an integrator instance */
    double node_rate;
    n.param<double>("node_rate", node_rate, 200);

    double sim_speed;
    n.param<double>("simulation_speed", sim_speed, 1.0);

    /** cast to seconds and round to ms */
    double dt = (1.0 / node_rate);
    dt = std::roundf(sim_speed * dt * 1000) / 1000;

    Dict params({{"tf",     dt},
                 {"tol",    1e-6},
                 {"method", CVODES}});
    Function ode = kite.getNumericDynamics();
    auto dynamics = kite.getAeroDynamicForces();

    ODESolver odeSolver(ode, params);

    Simulator simulator(odeSolver, n);
    simulator.sim_tether = simulate_tether;
    simulator.setNumericAirspeedMeas(kite.getNumericAirspeedMeas());
    simulator.setNumericAeroValues(kite.getNumericAeroValues());
    simulator.setNumericSpecNongravForce(kite.getNumericSpecNongravForce());
    simulator.setNumericSpecTethForce(kite.getNumericSpecTethForce());

    ros::Rate loop_rate(node_rate);

    while (ros::ok()) {
        if (simulator.is_initialized()) {
            simulator.simulate();

//            double sim_time_sec = (ros::Time::now() - simulation_time_start).toSec() * sim_speed;
            ros::Time sim_time;
//            sim_time.fromSec(sim_time_sec);
            sim_time = ros::Time::now();

            if (broadcast_state)
                simulator.publish_state(sim_time);
            else
                simulator.publish_pose(sim_time);

            if (-simulator.getState().nonzeros()[8] < 0) {
                std::cout << "\n----------------------------------\n"
                             "Hit the ground. End of simulation.\n"
                          << "----------------------------------\n";
                return 0;
            }
        }

        ros::spinOnce();
        loop_rate.sleep();
    }

    return 0;
}
