#include "casadi/casadi.hpp"
#include "iostream"

using namespace casadi;

int main(void)
{
    /** create basis functions */
    SX fx = 0;
    SX fy = 0;
    SX t = SX::sym("t");
    int K = 10;
    SX ax = SX::sym("ax", K);
    SX bx = SX::sym("bx", K);
    SX ay = SX::sym("ay", K);
    SX by = SX::sym("by", K);

    for(int k = 0; k < K; ++k)
    {
        fx += ax(k) * cos(k * t) + bx(k) * sin(k * t);
        fy += ay(k) * cos(k * t) + by(k) * sin(k * t);
    }

    SX f = SX::vertcat({fx,fy});
    SX p = SX::vertcat({ax,bx,ay,by});     // optimization variable
    Function F = Function("F",{p,t},{f});  // manifold

    double r = 10;
    SX fr = SX::vertcat({r * cos(t), r * sin(t)});
    Function Fr = Function("Fr",{t},{fr});

    /** generate obstacle samples */
    DMVector obstacles;
    double dx = 0.2;
    double dy = 0.1;
    // first obstacle
    double x = -7;
    double y = -1;
    while(x <= -2.5)
    {
        DM tmp = DM::vertcat({x,-1.0});
        obstacles.push_back(tmp);
        x += dx;
    }
    while(y <= 1)
    {
        DM tmp = DM::vertcat({-2.5,y});
        obstacles.push_back(tmp);
        y += dy;
    }

    //second obstacle
    x = 2.5;
    y= -1;

    while(x <= 7)
    {
        DM tmp = DM::vertcat({x,-1.0});
        obstacles.push_back(tmp);
        x += dx;
    }
    while(y <= 1)
    {
        DM tmp = DM::vertcat({2.5,y});
        obstacles.push_back(tmp);
        y += dy;
    }

    std::cout << obstacles.size() << "\n";

    /** create an optimization problem */
    SX opt_var = p;
    SX cost = 0;
    SX align_constr, inclusion_constr;

    double gamma = 0.01;

    for(int i = 0; i < obstacles.size(); ++i)
    {
        // compute cost and residual
        SX ti = atan2(obstacles.at(i)(1), obstacles.at(i)(0));
        SX fr_ti = Fr({ti})[0];
        SX f_ti  = F(SXVector{p,ti})[0];
        SX residual = f_ti - fr_ti;
        cost += SX::dot(residual, residual);

        // constraints: ray alignment
        SX Lti = f_ti + fr_ti;
        SX align = SX::dot(Lti, fr_ti);
        align_constr = SX::vertcat({align_constr, align});

        // constraint : inclusion
        SX pi = obstacles.at(i);
        SX inclusion = SX::dot(pi - Lti, fr_ti);
        inclusion_constr = SX::vertcat({inclusion_constr, inclusion});
    }
    cost += gamma * SX::dot(p,p);

    /** set the bounds */
    DM lbx = -DM::inf(opt_var.size1());
    DM ubx = DM::inf(opt_var.size1());

    DM lb_align = DM::zeros(align_constr.size1());
    DM ub_align = lb_align;

    DM lb_incl = DM::zeros(inclusion_constr.size1());
    DM ub_incl = DM::inf(inclusion_constr.size1());

    SX constraints = SX::vertcat({align_constr, inclusion_constr});
    DM lbg = DM::vertcat({lb_align, lb_incl});
    DM ubg = DM::vertcat({ub_align, ub_incl});

    /** create an NLP */
    SXDict QP;
    Function QP_Solver;
    Dict OPTS;
    DMDict ARG;

    QP["x"] = opt_var;
    QP["f"] = cost;
    QP["g"] = constraints;

    /** solver options */
    OPTS["ipopt.linear_solver"]         = "ma97";
    OPTS["ipopt.print_level"]           = 5;
    OPTS["ipopt.tol"]                   = 1e-4;
    OPTS["ipopt.acceptable_tol"]        = 1e-4;
    OPTS["ipopt.max_iter"]              = 5000;
    //OPTS["ipopt.warm_start_init_point"] = "yes";
    //OPTS["ipopt.hessian_approximation"] = "limited-memory";

    QP_Solver = casadi::nlpsol("solver", "ipopt", QP, OPTS);

    /** set default args */
    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;

    /** solve the problem */
    DMDict res = QP_Solver(ARG);
    DM qp_x    = res.at("x");

    std::cout << qp_x << "\n";

    return 0;
}
