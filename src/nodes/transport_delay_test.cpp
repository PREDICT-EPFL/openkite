#include "ros/ros.h"
#include "test_package/aircraft_controls.h"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/thread.hpp"
#include <random>

/** PRBS generator node */
int main(int argc, char** argv)
{
    ros::init(argc, argv, "transport_delay_test");
    ros::NodeHandle handle;
    ros::Rate loop_rate(30);

    ros::Publisher pub = handle.advertise<test_package::aircraft_controls>("/controls", 100);
    test_package::aircraft_controls msg;
    msg.ailerons = 0;
    msg.flaps    = 0;
    msg.rudder   = 1;
    msg.thrust   = 0;

    double rudder = 1;

    /** PRBS generator initialization */
    std::random_device device;
    std::mt19937       generator(device());
    std::uniform_int_distribution<> distribution(20, 50);

    int counter = 0;

    while(ros::ok())
    {
        //        int duty_cycle = distribution(generator);
        //        rudder = (rudder == 1) ? 0 : 1;
        //        msg.rudder = rudder;
        //        msg.header.stamp = ros::Time::now();
        //        pub.publish(msg);
        //        boost::this_thread::sleep(boost::posix_time::milliseconds(duty_cycle));

        if (counter < 200)
        {
            msg.rudder = 1.0;
            //std::cout << "case1: " << counter << "\n";
        }
        else if ((counter >= 200) && (counter < 400))
        {
            msg.rudder = 0.0;
            //std::cout << "case2: " << counter << "\n";
        }
        else
            counter = 0;

        msg.header.stamp = ros::Time::now();
        pub.publish(msg);
        ++counter;

        loop_rate.sleep();
    }

    return 0;
}
