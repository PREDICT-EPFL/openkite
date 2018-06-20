#ifndef BOAT_MATH_H
#define BOAT_MATH_H

#include <casadi/casadi.hpp>

namespace boat_math {

casadi::SX quatmul(const casadi::SX &q1, const casadi::SX &q2);
casadi::SX quatconj(const casadi::SX &q);
casadi::SX quatinv(const casadi::SX &q);
casadi::SX quatrot(const casadi::SX &q, const casadi::SX &r);
casadi::SX quatrot_inverse(const casadi::SX &q, const casadi::SX &r);
casadi::SX rk4_symbolic(const casadi::SX &x,
                        const casadi::SX &u,
                        casadi::Function &fuc,
                        const casadi::SX &h);
casadi::SX skew_mat(const casadi::SX &v);

}

#endif /* BOAT_MATH_H */
