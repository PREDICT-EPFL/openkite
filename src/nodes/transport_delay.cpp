#include "ros/ros.h"
#include "boost/thread.hpp"
#include "boost/thread/mutex.hpp"
#include "boost/circular_buffer.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "openkite/aircraft_controls.h"
#include <random>

class TDelay
{
public:
    TDelay(const ros::NodeHandle &n);
    virtual ~TDelay(){}

    /** main loop */
    int  getNumPublishers(){return sub.getNumPublishers();}
    void run_publisher();

private:
    double delay_mean;
    double delay_deviation;
    double loop_rate;

    std::random_device rndev;
    std::mt19937 generator;
    std::uniform_real_distribution<> distribution;

    boost::circular_buffer<openkite::aircraft_controls> buffer;
    boost::thread pub_thread;
    boost::mutex  m_mutex;

    std::shared_ptr<ros::NodeHandle> nh;
    ros::Publisher  pub;
    ros::Subscriber sub;

    void inputCallback(const openkite::aircraft_controls::Ptr &input);
    void publish();
};

TDelay::TDelay(const ros::NodeHandle &n)
{
    nh = std::make_shared<ros::NodeHandle>(n);

    nh->param<double>("delay", delay_mean, 20.0);
    nh->param<double>("deviation", delay_deviation, 5.0);
    nh->param<double>("rate", loop_rate, 5.0);

    pub = nh->advertise<openkite::aircraft_controls>("/delayed_control", 100);
    sub = nh->subscribe("/controls", 100, &TDelay::inputCallback, this);

    pub_thread = boost::thread(boost::bind(&TDelay::publish, this));
    buffer.set_capacity(2);

    generator = std::mt19937(rndev());
    double lower_bound = (delay_mean - delay_deviation > 0) ? (delay_mean - delay_deviation) : 0;
    double upper_bound = delay_mean + delay_deviation;
    distribution     = std::uniform_real_distribution<>(lower_bound, upper_bound);
}

void TDelay::inputCallback(const openkite::aircraft_controls::Ptr &input)
{
    boost::mutex::scoped_try_lock lock(m_mutex);
    if(lock)
    {
        buffer.push_back(*input);
    }
}

void TDelay::run_publisher()
{
    pub_thread.join();
}

void TDelay::publish()
{
    openkite::aircraft_controls msg;
    while (ros::ok())
    {
        if(!buffer.empty())
        {
            double _delay = distribution(generator);
            {
                /** take an element of the queue publish it and delete */
                boost::mutex::scoped_lock lock(m_mutex);
                msg = buffer.front();
                buffer.pop_front();
            }
            boost::this_thread::sleep(boost::posix_time::milliseconds(_delay));
            msg.header.stamp = ros::Time::now();
            pub.publish(msg);
        }
        else
        {
            boost::this_thread::sleep(boost::posix_time::milliseconds(1000 / loop_rate ));
        }
    }
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "transport_delay");
    ros::NodeHandle handle;

    double rate;
    handle.param<double>("rate", rate, 40);

    TDelay transport_delay(handle);
    ros::Rate loop_rate(rate);

    while(ros::ok())
    {
        if(transport_delay.getNumPublishers())
        {
            ros::spinOnce();
            loop_rate.sleep();
        }
        else
        {
            loop_rate.sleep();
        }
    }

    transport_delay.run_publisher();

    return 0;
}
