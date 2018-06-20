#include <casadi/casadi.hpp>
using namespace casadi;

namespace boat_math {

SX quatmul(const SX &q1, const SX &q2)
{
    SX s1 = q1(0);
    SX v1 = q1(Slice(1,4));

    SX s2 = q2(0);
    SX v2 = q2(Slice(1,4));

    SX q = SX::vertcat({s1*s2 - SX::dot(v1,v2), SX::cross(v1,v2) + s1 * v2 + s2 * v1});
    return q;
}

SX quatconj(const SX &q)
{
    return SX::vertcat({q(0), -q(1), -q(2), -q(3)});
}

SX quatinv(const SX &q)
{
    // TODO: should be divided by squared norm
    return quatconj(q) / SX::norm_2(q);

    // SX norm_squared = SX::dot(q,q);
    // return quatconj(q) / norm_squared;
}

// TODO: verify correctness (multiplication order)
SX quatrot(const SX &q, const SX &r)
{
    // SX qinv = quatinv(q);
    SX qinv = quatconj(q);  // q_inv = q_conj for unit quaternions
    SX qr   = SX::vertcat({0, r});
    SX qrr  = quatmul(quatmul(qinv, qr),q);
    return qrr(Slice(1,4));
}

SX quatrot_inverse(const SX &q, const SX &r)
{
    // SX qinv = quatinv(q);
    SX qinv = quatconj(q);  // q_inv = q_conj for unit quaternions
    SX qr   = SX::vertcat({0, r});
    SX qrr  = quatmul(quatmul(q, qr),qinv);
    return qrr(Slice(1,4));
}

SX rk4_symbolic(const SX &x,
                const SX &u,
                Function &func,
                const SX &h)
{
    SXVector res = func(SXVector{x, u});
    SX k1 = res[0];
    res = func(SXVector{x + 0.5 * h * k1, u});
    SX k2 = res[0];
    res = func(SXVector{x + 0.5 * h * k2, u});
    SX k3 = res[0];
    res = func(SXVector{x + h * k3, u});
    SX k4 = res[0];

    return x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
}

SX skew_mat(const SX &v)
{
    return SX::vertcat({
        SX::horzcat({0, -v(2), v(1)}),
        SX::horzcat({v(2), 0, -v(0)}),
        SX::horzcat({-v(1), v(0), 0})
    });
}

} // namespace boat_math