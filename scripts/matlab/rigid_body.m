function [model] = rigid_body(dt)
%rigid body kinematic equations
import casadi.*

r = SX.sym('r', 3);
q = SX.sym('q', 4);
vb = SX.sym('vb',3);
wb = SX.sym('wb',3);
u = SX.sym('u',3);

%translation
qv = quat_mul(q, [0; vb]);
qv_q = quat_mul(qv, quat_inv(q));
rdot = qv_q(2:4);
%rotation
lambda = -10;
qdot = 0.5 * quat_mul(q, [0; wb]) + 0.5 * lambda * q * (q' * q - 1);

%velocities
vb_dot = zeros(3,1);
wb_dot = zeros(3,1);

%define reduced integrator
kinematics_red = [rdot; qdot];
state_red = [r; q];
params_red = [vb; wb];

ode_red = struct('x',state_red, 'p', params_red, 'ode', kinematics_red);
opts_red = struct('tf', dt);
INT_RED = casadi.integrator('CVODES_INT','cvodes', ode_red, opts_red);
model.INT_RED = INT_RED;

%full model
state = [vb; wb; r; q];
kinematics = [vb_dot; wb_dot; rdot; qdot];
%jacobian
jacobian = kinematics.jacobian(state);
JACOBIAN = Function('RB_JACOBIAN',{state, u},{jacobian});
model.JACOBIAN = JACOBIAN;

%integration
ode = struct('x',state, 'ode',kinematics);
opts = struct('tf', dt);
INTEGRAL = casadi.integrator('CVODES_INT','cvodes', ode, opts);
model.INT = INTEGRAL;
end


function [q] = quat_mul(q1,q2)
s1 = q1(1);
v1 = q1(2:4);

s2 = q2(1);
v2 = q2(2:4);

s = s1*s2 - dot(v1,v2);
v = cross(v1,v2) + s1 * v2 + s2 * v1;

q = [s; v(1); v(2); v(3)];
end

%inverse quaternion
function [q_inv] = quat_inv(q)
q_inv = [q(1); -q(2); -q(3); -q(4)];
end

%L2-norm of a vector
function [L2] = L2(x)
if(size(x,1) >= size(x,2))
    L2 = sqrt(x' * x);
else
    L2 = sqrt(x * x');
end

end