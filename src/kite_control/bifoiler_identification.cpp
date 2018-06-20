#include "kiteNMPF.h"
#include "integrator.h"

#define BOOST_TEST_MODULE bifoiler_identification
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include "pseudospectral/chebyshev.hpp"

#include "boat_model.h"

using namespace casadi;

BOOST_AUTO_TEST_SUITE( bifoiler_identification )

BOOST_AUTO_TEST_CASE( bifoiler_id_test )
{
    /** define boat dynamics */
    std::string boat_config_file = "config.yaml";
    BoatProperties boat_props = BoatProperties::Load(boat_config_file);

    bifoiler::BoatDynamics boat(boat_props);
    bifoiler::BoatDynamics boat_int(boat_props); //integration model
    Function ode = boat_int.getNumericDynamics();

    /** get dynamics function and state Jacobian */
    Function DynamicsFunc = boat.getNumericDynamics();
    SX X = boat.getSymbolicState();
    SX U = boat.getSymbolicControl();
    //SX P = boat.getSymbolicParameters(); not there so far

    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_SUITE_END()
