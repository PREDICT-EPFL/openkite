#include "homotopy/psarc.hpp"

struct FX
{
    FX();
    ~FX(){}

    using num = casadi::DM;
    using sym = casadi::SX;

    /** symbolic variable */
    sym var;
    sym sym_func;
    sym sym_jac;

    casadi::Function func;
    casadi::Function jac;

    num eval(const num &arg)
    {
        return func({arg})[0];
    }

    num jacobian(const num &arg)
    {
        return jac({arg})[0];
    }

    sym operator()()
    {
        return sym_func;
    }
};

FX::FX()
{
    var = casadi::SX::sym("x");
    sym_func = var + 4;
    sym_jac = casadi::SX::jacobian(sym_func, var);

    func = casadi::Function("lox",{var},{sym_func});
    jac = casadi::Function("pidor",{var},{sym_jac});
}

int main(void)
{	
    casadi::DM x0 = 1;
    symbolic_psarc<FX, casadi::Dict>psarc(x0);
    return 0;
}
