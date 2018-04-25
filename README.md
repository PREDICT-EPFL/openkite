# openKITE
OpenKITE is a ROS package for simulation, estimation and control of rigid-wing airborne wind energy kites. Developed at the
Automatic Control laboratory EPFL as part of the AWESCO project.

-- **Dependencies**:
 - CasADi (v3.0.0-rc2) : https://github.com/casadi
 - YAMLCPP (v0.5.3) : https://github.com/jbeder/yaml-cpp
 - ROS (Jade, Kinetic, Lunar) : http://www.ros.org/
 
-- **Build**:

OpenKITE is based on the Catkin building system (http://wiki.ros.org/catkin). In order to build the package place it in your 
Catkin workspace and call:
```
cd your_catkin_ws
catkin_make
```

-- **Content**:
 - *kite_model* : contains collection of kite models, integrators and flight simulator code;
 - *kite_control* : contains implementation of a Nonlinear Model Predictive controller for orbit tracking, and kite 
 identification tool;
 - *kite_estimation* : contains imlpementation of the Extended Kalman filter for state estimation;
 - *scripts* : collection of Matlab and Python scripts for log/data analysis and calibration;
 - *launch* : ROS launch files to ease using the software;
 - *data* : sample kite configuration file;
 - *arduino* : transmit control commands over radio link (Spectrum receiver/transmitter) using ros-serial framework
 - **doc** : user manual
 
 -- **Contact**: plistov@gmail.com
