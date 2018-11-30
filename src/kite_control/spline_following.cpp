
#include "kiteNMPF.h"

using namespace casadi;

SX polynomial(const SX &coef, const SX &x0, const SX &x)
{
    SX poly = SX::vertcat(SXVector{pow(x-x0, 3), pow(x-x0, 2), x-x0, 1});
    return SX::dot(coef, poly);
}


int main(void)
{
    SX breaks = DM({0, 1.4225, 2.8451, 4.2676, 5.6902, 7.1127, 8.5353, 9.9578, 11.3803});
    SX coeffs_x = DM::horzcat(DMVector{DM({-0.0049,    0.0510,   -0.0261,   -5.0000}),
                                     DM({-0.0049,    0.0300,    0.0890,   -4.9482}),
                                     DM({0.0555,    0.0090,    0.1445,   -4.7750}),
                                     DM({-0.0635,    0.2459,    0.5071,   -4.3915}),
                                     DM({-0.0173,   -0.0252,    0.8211,   -3.3554}),
                                     DM({0.0212,   -0.0989,    0.6446,   -2.2881}),
                                     DM({-0.0003,   -0.0082,    0.4923,   -1.5100}),
                                     DM({-0.0003,   -0.0093,    0.4674,   -0.8270}) });

    SX coeffs_y = DM::horzcat(DMVector{DM({0.0005,   -0.0047,    1.0048,   -0.0000}),
                                       DM({0.0005,   -0.0025,    0.9945,    1.4213}),
                                       DM({-0.0148,   -0.0002,    0.9908,    2.8326}),
                                       DM({-0.1311,   -0.0632,    0.9007,    4.1992}),
                                       DM({0.1598,   -0.6226,   -0.0749,    4.9753}),
                                       DM({-0.0221,    0.0594,   -0.8761,    4.0689}),
                                       DM({0.0071,   -0.0351,   -0.8415,    2.8790}),
                                       DM({0.0071,   -0.0050,   -0.8984,    1.6313}) });

    SX x = SX::sym("x");
    coeffs_x = coeffs_x.T();
    coeffs_y = coeffs_y.T();

    SX path_x = 0;
    SX path_y = 0;
    for(int i = 0; i <= 7; ++i)
    {
        SX pix = polynomial(coeffs_x(i, Slice(0, 4)).T(), SX(breaks[i]), x);
        path_x += pix * (heaviside(x - breaks[i] + 1e-5) - heaviside(x - breaks[i+1])); // 1e-5 - Heaviside correction

        SX piy = polynomial(coeffs_y(i, Slice(0, 4)).T(), SX(breaks[i]), x);
        path_y += piy * (heaviside(x - breaks[i] + 1e-5) - heaviside(x - breaks[i+1])); // 1e-5 - Heaviside correction
    }

    Function path_funx = Function("path", {x}, {path_x});
    Function path_funy = Function("path", {x}, {path_y});
    std::vector<double> vecx;
    std::vector<double> vecy;
    std::vector<double> vect;
    double length = 0.0;
    while(length <= 11.0)
    {
        /// x
        DM res = path_funx(DMVector{length})[0];
        vecx.push_back(res.nonzeros()[0]);
        /// y
        res = path_funy(DMVector{length})[0];
        vecy.push_back(res.nonzeros()[0]);

        vect.push_back(length);
        length += 0.1;
    }

    std::cout << vect << "\n";
    std::cout << vecx << "\n";
    std::cout << vecy << "\n";

    return 0;
}
