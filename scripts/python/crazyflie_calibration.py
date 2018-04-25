#!/usr/bin/env python
import roslib
import rospy
import numpy as np

from geometry_msgs.msg import Twist
from test_package.msg import aircraft_controls
from threading import Thread

import sys, select, termios, tty

thrust_inc = 200
angle_inc  = 200

thrust   = 0
elevator = 0
rudder   = 0
status = 0


speedBindings={
		'q':(thrust_inc,0,0),
		'z':(-thrust_inc,0,0),
		'w':(0,angle_inc,0),
		'x':(0,-angle_inc,0),
		'e':(0,0,angle_inc),
		'c':(0,0,-angle_inc),
	      }

def getKey():
    tty.setraw(sys.stdin.fileno())
    select.select([sys.stdin], [], [], 0)
    key = sys.stdin.read(1)
    termios.tcsetattr(sys.stdin, termios.TCSADRAIN, settings)
    return key

def process_keyboard():
    try:
        print vels(0,0,0)
        while(1):
            key = getKey()
            if key in speedBindings.keys():
                global thrust, elevator, rudder, status
                thrust = thrust + speedBindings[key][0]
                elevator = elevator + speedBindings[key][1]
                rudder  = rudder + speedBindings[key][2]

                print vels(thrust,elevator,rudder)
                status = (status + 1) % 15
            else:
                if (key == '\x03'):
                    break

            #twist = Twist()
            #twist.linear.x = elevator
            #twist.linear.y = rudder
            #twist.linear.z = thrust


    except Exception, e:
        print e

    finally:
        twist = aircraft_controls()
        twist.thrust = 0
        twist.rudder = 0
        twist.elevator = 0
        pub.publish(twist)
        termios.tcsetattr(sys.stdin, termios.TCSADRAIN, settings)



def vels(thrust, elevator, rudder):
	return "currently:\tthrust %s\televator %s\trudder %s" % (thrust,elevator,rudder)

if __name__=="__main__":
    settings = termios.tcgetattr(sys.stdin)

    #pub = rospy.Publisher('/crazyflie/cmd_vel', Twist, queue_size = 1)
    pub = rospy.Publisher('/kite_controls', aircraft_controls, queue_size = 1)
    rospy.init_node('talker')
    rate = rospy.Rate(20)

    t = Thread(target=process_keyboard)
    t.start()

    while not rospy.is_shutdown():
        twist = aircraft_controls()
        twist.thrust = thrust
        twist.rudder = rudder
        twist.elevator = elevator

        #twist.angular.x = 0; twist.angular.y = 0; twist.angular.z = 0
        pub.publish(twist)
        rate.sleep()



