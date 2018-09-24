#ifndef PSARC_HPP
#define PSARC_HPP

#include "casadi/casadi.hpp"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/Eigenvalues"

namespace psarc_math
{

enum MAT : int {CASADI, EIGEN_DENSE, EIGEN_SPARSE};

static casadi::DM solve(const casadi::DM &A, const casadi::DM &b, MAT mat_type = MAT::EIGEN_DENSE)
{
    casadi::DM x;
    /** if A dimansion is small use native Casadi solver */
    switch(mat_type){
    case CASADI : {
        x = casadi::DM::solve(A, b);
        return x;
    }
    case EIGEN_DENSE : {
        Eigen::MatrixXd _A;
        Eigen::MatrixXd _b;
        _A = Eigen::MatrixXd::Map(casadi::DM::densify(A).nonzeros().data(), A.size1(), A.size2());
        _b = Eigen::VectorXd::Map(casadi::DM::densify(b).nonzeros().data(), b.size1());

        /** solve the linear system and cast back to Casadi types */
        Eigen::VectorXd sol = _A.partialPivLu().solve(_b);
        std::vector<double> sold;
        sold.resize(static_cast<size_t>(sol.size()));
        Eigen::Map<Eigen::VectorXd>(&sold[0], sol.size()) = sol;
        x = casadi::DM(sold);
        return x;
    }
    case EIGEN_SPARSE : {
        Eigen::SparseMatrix<double> _A = C2ESparse<double>(A);
        return casadi::DM();
    }
    }

}

template<typename Scalar>
Eigen::SparseMatrix<Scalar> C2ESparse(const casadi::DM &matrix)
{
    casadi::Sparsity SpA = matrix.get_sparsity();
    std::vector<int> output_row, output_col;
    SpA.get_triplet(output_row, output_col);
    std::vector<double> values = matrix.get_nonzeros();

    using T = Eigen::Triplet<double>;
    std::vector<T> TripletList;
    TripletList.resize(values.size());
    for(int k = 0; k < values.size(); ++k)
        TripletList[k] = T(output_row[k], output_col[k], values[k]);

    Eigen::SparseMatrix<double> SpMatrx(matrix.size1(), matrix.size2());
    SpMatrx.setFromTriplets(TripletList.begin(), TripletList.end());

    return SpMatrx;
}

// end of namespace
}


template<typename FX, typename CorrectorProps>
typename FX::x psarc(const typename FX::x &init_guess, const CorrectorProps &props = CorrectorProps())
{
    FX system;

    typename FX::x res = system.eval(init_guess);
    std::cout << res << "\n";

    res = system.jac(init_guess);
    std::cout << res << "\n";
    return res;
}

template<typename Equalities, typename CorrectorProps>
class symbolic_psarc
{
public:
    symbolic_psarc(const typename Equalities::num &init_guess, const CorrectorProps &props = CorrectorProps());
    ~symbolic_psarc(){}

    casadi::DMDict operator()();
};

template<typename Equalities, typename CorrectorProps>
symbolic_psarc<Equalities, CorrectorProps>::symbolic_psarc(const typename Equalities::num &init_guess, const CorrectorProps &props)
{
    /** create an instance of the system */
    Equalities FX;

    /** generate convex homotopy equation */
    typename Equalities::sym x = FX.var;
    typename Equalities::sym lambda = Equalities::sym::sym("lambda");
    typename Equalities::sym x0 = {init_guess};
    typename Equalities::sym homotopy = (lambda) * (x - x0) + (1 - lambda) * FX();

    /** create a corrector */
    x = Equalities::sym::vertcat({x, lambda});
    typename Equalities::sym w = Equalities::sym::sym("w", x.size1(), 1);

    casadi::SXDict NLP;
    casadi::Dict   OPTS;
    casadi::DMDict ARG;
    NLP["x"] = x;
    NLP["f"] = 0.5 * Equalities::sym::dot((x - w), (x - w));
    NLP["g"] = homotopy;
    NLP["p"] = w;

    OPTS["ipopt.linear_solver"]  = "ma97";
    OPTS["ipopt.print_level"]    = 1;
    OPTS["ipopt.tol"]            = 1e-6;
    OPTS["ipopt.acceptable_tol"] = 1e-6;
    OPTS["ipopt.warm_start_init_point"] = "yes";

    casadi::Function Corrector = casadi::nlpsol("solver", "ipopt", NLP, OPTS);

    typename Equalities::num lbx = -Equalities::num::inf(x.size1());
    typename Equalities::num ubx = Equalities::num::inf(x.size1());
    typename Equalities::num lbg = Equalities::num::zeros(homotopy.size1());
    typename Equalities::num ubg = lbg;

    ARG["lbx"] = lbx;
    ARG["ubx"] = ubx;
    ARG["lbg"] = lbg;
    ARG["ubg"] = ubg;
    ARG["p"] = Equalities::num::vertcat({1,1});

    ARG["lbx"](lbx.size1() - 1, 0) = 1;
    ARG["ubx"](ubx.size1() - 1, 0) = 1;

    /** solve initial homotopy equations */
    casadi::DMDict solution = Corrector(ARG);
    casadi::Dict stats = Corrector.stats();

    std::cout << solution.at("x") << "\n";

    /** implement predict-corrector scheme */
    typename Equalities::sym homotopy_jac = Equalities::sym::jacobian(homotopy, x);
    casadi::Function h_jac_func = casadi::Function("h_jac", {x}, {homotopy_jac});
    double lambda_val = 1.0;
    double h = 0.02;
    typename Equalities::num x_next;
    typename Equalities::num x0_num = solution.at("x");

    ARG["lbx"](lbx.size1() - 1, 0) = -Equalities::num::inf(1);
    ARG["ubx"](ubx.size1() - 1, 0) = Equalities::num::inf(1);;

    while(lambda_val > 0.0)
    {
        /** estimate tangent direction*/
        /** ... */

        /** apply corrector */
        ARG["p"] = x_next;
        //ARG["lbx"](lbx.size1() - 1, 0) = x_next(1);
        //ARG["ubx"](ubx.size1() - 1, 0) = x_next(1);
        solution = Corrector(ARG);

        /** prepare for the next step */
        x0_num = solution.at("x");

        if(x0_num(1).nonzeros()[0] < 0.0)
        {   /** refine solution if it crosses 0 */
            x0_num(1) = 0.0;
            ARG["lbx"](lbx.size1() - 1, 0) = x0_num(1);
            ARG["ubx"](ubx.size1() - 1, 0) = x0_num(1);
            solution = Corrector(ARG);
            x0_num = solution.at("x");
        }
        std::cout << "Corrector: " << x0_num << "\n";
        std::cout << Corrector.stats() << "\n";

        lambda_val = x0_num(1).nonzeros()[0];
    }
}

#endif // PSARC_HPP
