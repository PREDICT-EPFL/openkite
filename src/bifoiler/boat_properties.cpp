#include <yaml-cpp/yaml.h>

#include "boat_properties.h"

BoatProperties BoatProperties::Load(const std::string &filename)
{
    YAML::Node config = YAML::LoadFile(filename);

    BoatProperties prop;
    prop.name = config["name"].as<std::string>();

    auto estimator = config["estimator"];
    prop.estimator.t_samp = estimator["t_samp"].as<double>();

    auto env = config["env"];
    prop.env.rho_sh2o = env["rho_sh2o"].as<double>();
    prop.env.g = env["g"].as<double>();

    auto sensor = config["sensor"];
    prop.sensor.r_ant[0] = sensor["r_ant"][0].as<double>();
    prop.sensor.r_ant[1] = sensor["r_ant"][1].as<double>();
    prop.sensor.r_ant[2] = sensor["r_ant"][2].as<double>();

    auto inertia = config["inertia"];
    prop.inertia.mass = inertia["mass"].as<double>();
    prop.inertia.mass_cad = inertia["mass_cad"].as<double>();
    prop.inertia.Ixy = inertia["Ixy"].as<double>();
    prop.inertia.Ixz = inertia["Ixz"].as<double>();
    prop.inertia.Iyz = inertia["Iyz"].as<double>();
    prop.inertia.Ixx = inertia["Ixx"].as<double>();
    prop.inertia.Iyy = inertia["Iyy"].as<double>();
    prop.inertia.Izz = inertia["Izz"].as<double>();

    auto foils = config["foils"];
    prop.foils.ARff = foils["ARff"].as<double>();
    prop.foils.areaff = foils["areaff"].as<double>();
    prop.foils.mac = foils["mac"].as<double>();
    prop.foils.wingspanff = foils["wingspanff"].as<double>();

    auto hydro = config["hydrodynamic"];
    prop.hydrodynamic.CL0 = hydro["CL0"].as<double>();
    prop.hydrodynamic.CLa_total = hydro["CLa_total"].as<double>();
    prop.hydrodynamic.e_oswald = hydro["e_oswald"].as<double>();
    prop.hydrodynamic.CD0_total = hydro["CD0_total"].as<double>();
    prop.hydrodynamic.CYb = hydro["CYb"].as<double>();
    prop.hydrodynamic.Cm0 = hydro["Cm0"].as<double>();
    prop.hydrodynamic.Cma = hydro["Cma"].as<double>();
    prop.hydrodynamic.Cn0 = hydro["Cn0"].as<double>();
    prop.hydrodynamic.Cnb = hydro["Cnb"].as<double>();
    prop.hydrodynamic.Cl0 = hydro["Cl0"].as<double>();
    prop.hydrodynamic.Clb = hydro["Clb"].as<double>();
    prop.hydrodynamic.CLq = hydro["CLq"].as<double>();
    prop.hydrodynamic.Cmq = hydro["Cmq"].as<double>();
    prop.hydrodynamic.CYr = hydro["CYr"].as<double>();
    prop.hydrodynamic.Cnr = hydro["Cnr"].as<double>();
    prop.hydrodynamic.Clr = hydro["Clr"].as<double>();
    prop.hydrodynamic.CYp = hydro["CYp"].as<double>();
    prop.hydrodynamic.Clp = hydro["Clp"].as<double>();
    prop.hydrodynamic.Cnp = hydro["Cnp"].as<double>();
    prop.hydrodynamic.CXdf = hydro["CXdf"].as<double>();
    prop.hydrodynamic.CYdr = hydro["CYdr"].as<double>();
    prop.hydrodynamic.CZde = hydro["CZde"].as<double>();
    prop.hydrodynamic.CZdf = hydro["CZdf"].as<double>();
    prop.hydrodynamic.CLda = hydro["CLda"].as<double>();
    prop.hydrodynamic.CLdr = hydro["CLdr"].as<double>();
    prop.hydrodynamic.CMde = hydro["CMde"].as<double>();
    prop.hydrodynamic.CMdf = hydro["CMdf"].as<double>();
    prop.hydrodynamic.CNda = hydro["CNda"].as<double>();
    prop.hydrodynamic.CNdr = hydro["CNdr"].as<double>();

    return prop;
}
