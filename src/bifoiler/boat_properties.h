#ifndef BOAT_PROPERTIES_H
#define BOAT_PROPERTIES_H

struct BoatProperties
{
    std::string name;

    struct {
        double t_samp;
    } estimator;

    struct {
        double rho_sh2o;
        double g;
    } env;

    struct {
        double r_ant[3];
    } sensor;

    struct {
        double mass;
        double mass_cad;
        double Ixy;
        double Ixz;
        double Iyz;
        double Ixx;
        double Iyy;
        double Izz;
    } inertia;

    struct {
        double ARff;
        double areaff;
        double mac;
        double wingspanff;
    } foils;

    struct {
        double CL0;
        double CLa_total;
        double e_oswald;
        double CD0_total;
        double CYb;
        double Cm0;
        double Cma;
        double Cn0;
        double Cnb;
        double Cl0;
        double Clb;
        double CLq;
        double Cmq;
        double CYr;
        double Cnr;
        double Clr;
        double CYp;
        double Clp;
        double Cnp;
        double CXdf;
        double CYdr;
        double CZde;
        double CZdf;
        double CLda;
        double CLdr;
        double CMde;
        double CMdf;
        double CNda;
        double CNdr;
    } hydrodynamic;

    static BoatProperties Load(const std::string &filename);
    // static void Save(const std::string &filename); // TODO
};

#endif /* BOAT_PROPERTIES_H */
