#include <control_proxy_node.hpp>

void ControlProxyNode::controlCallback(const openkite::aircraft_controls::ConstPtr& msg)
{
    /** transform data */
    servo_msg.data[0] = static_cast<int16_t> (1100 + (800 / 0.15) * msg->thrust);
    servo_msg.data[1] = static_cast<int16_t> (1500 + (400 / 0.26) * msg->rudder);
    servo_msg.data[2] = static_cast<int16_t> (1500 + (400 / 0.26) * msg->elevator);
    servo_msg.data[3] = static_cast<int16_t> (1500 + (400 / 0.26) * msg->ailerons);
}

ControlProxyNode::ControlProxyNode(const ros::NodeHandle &_nh)
{
    nh = std::make_shared<ros::NodeHandle>(_nh);
    proxy_sub = nh->subscribe("kite_controls", 100, &ControlProxyNode::controlCallback, this);
    proxy_pub = nh->advertise<std_msgs::Int16MultiArray>("servo_controls", 100);

    /** array initialization */
    servo_msg.layout.dim.push_back(std_msgs::MultiArrayDimension());
    servo_msg.layout.dim[0].label = "lox";
    servo_msg.layout.dim[0].size = 4;
    servo_msg.layout.dim[0].stride = 4;
    servo_msg.layout.data_offset = 0;
    servo_msg.data.reserve(4);

    set_servos(1100, 1500, 1500, 1500);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "control_proxy_node");
    ros::NodeHandle n;
    ControlProxyNode proxy(n);

    ros::Rate loop_rate(50); /** 40 Hz */
    uint counter = 0;

    bool testing = false;

    while (ros::ok())
    {
        // testing
        if (testing)
        {
            if (counter < 200)
            {
                proxy.set_servos(1100, 1500, 1800, 1500);
                //std::cout << "case1: " << counter << "\n";
            }
            else if ((counter >= 200) && (counter < 400))
            {
                proxy.set_servos(1100, 1500, 1200, 1500);
                //std::cout << "case2: " << counter << "\n";
            }
            else
                counter = 0;

            ++counter;
        }

        proxy.publish();
        ros::spinOnce();
        loop_rate.sleep();
    }

    return 0;
}
