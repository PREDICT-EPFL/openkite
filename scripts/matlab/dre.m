function[DRE_INTEGRATOR] = dre(N, h)
%solve Differential Riccati Equation
% dX/dt = AX + XA' + Q
X = casadi.SX.sym('X', N,N);
A = casadi.SX.sym('A', N,N);
Q = casadi.SX.sym('Q', N,N);
%h = casadi.SX.sym('time_int', 1);

%create dynamics
Xdot = riccati(vec(X), A, Q);
%xdot_func = casadi.Function('riccati_ode', {vec(X), A, Q}, {Xdot});

%create integrator 
%integrator = RK4_DRE_SYM(vec(X), A, Q, h, xdot_func);
%DRE_INTEGRATOR = casadi.Function('dre_integrator', {vec(X), A, Q, h}, {integrator});
%CVODES integrator
params = [vec(A);vec(Q)];
ode = struct('x',vec(X), 'p',params, 'ode', Xdot);
opts = struct('tf', h);
DRE_INTEGRATOR = casadi.integrator('CVODES_INT','cvodes', ode, opts);
end

function[xdot] = riccati(x, A, Q)
X = reshape(x, size(A));
xdot = vec(A * X) + vec(X * A') + vec(Q);
end

%Runge-Kutta 4th-order symbolic integration method
function [X] = RK4_DRE_SYM(x, a, q, h, functor)
% this integrator is problem specific
k1 = functor({x, a, q});
k2 = functor({x + 0.5 * h * k1{1}, a, q});
k3 = functor({x + 0.5 * h * k2{1}, a, q});
k4 = functor({x + h * k3{1}, a, q});

x_h = x + (h/6) * (k1{1} + 2*k2{1} + 2*k3{1} + k4{1});

X = reshape(x_h, size(a));
end