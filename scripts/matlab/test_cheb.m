function [] = test_cheb()
%spectral integration test

%equation xdot - lambda * x = 0

clc;
close all;

N = 10;
[x,D] = cheb(N);
D(end,:) = zeros(N+1,1);
D(end, end) = 1;

u = [zeros(N, 1); 1];
lambda = -1.0;  % fixed-point iteration
I = eye(size(D));
I(end,end) = 0;

%solve
u = (D - lambda * I) \ u;

clf, subplot('position',[.1 .4 .8 .5])
plot(x,u,'.','markersize',16)
xx = 0:.01:1;
hold on
grid on
t = 0:0.1:1;
plot(t, exp(lambda * t), 'red')


%-------------------------------------------------------------------------%
%
% solve system of ODEs xdot = Ax + Bu
A = diag([-1;-2;-3]);
A(1,3) = -0.5;
B = [1;1;1];
C = eye(3,3);

%set up exact solution
options.time = 'range';
time = [0,1,0.01];
sys.a = A;
sys.b = B;
sys.c = C;
init = 10 * rand(3,1);

[X,Y] = ss_solver(sys, init, time, options);

%solve the system with Chebyshev collocation
S = diag([0.1;0.1;0.1]);
D3 = kron(D, eye(3));
Ma = kron(eye(N+1), S * A * inv(S));
Ma(end-2:end, end-2:end) = zeros(3);
pseudo_init = [init];
u = [repmat(S * B,N,1); S * pseudo_init];

%solve the system
tic
u = (D3 - Ma) \ u;
toc

u = kron(eye(N+1), inv(S)) * u;

figure
plot([0:0.01:1], X);
hold on
grid on

Xu = [];
for i = 0:N
    Xu = [Xu; u((3*i + 1):(3*i + 3))'];
end
plot(x, Xu, '--');
hold off

%-------------------------------------------------------------------------%
%solve nonlinear ODE xdot = f(x,u)
aircraft = ReadYaml('umx_radian.yaml');
parameters.simulation = 0;
parameters.plot = 0;
parameters.vr = 0;
parameters.int_type = 'cvodes';

disp(aircraft)
disp(parameters)

%get aircraft model
[num, log, sym] = flight_sim(aircraft, parameters);

%get dynamics function and state Jacobian
dynamics = num.DYN_FUNC;
dyn_jacobian = num.DYN_JACOBIAN;

x0 = [7;0;0;0;0;0;0;0;0;1;0;0;0];
u0 = [0;0;0];

%get integrator
INT = num.INT;

%reference integration loop Adams-Bashforth
time = 0;
sim_time = 1;
log = [x0'];
time_line = [time];
while time < sim_time    
    out = INT(struct('x0', x0, 'p', u0));
    x0 = full(out.xf);
    time = time + 0.02;
    time_line = [time_line; time];
    log = [log; x0'];
end

figure
plot(time_line, log(:,1:3));
grid on
hold on

%solve scaled system
%scaling matrix
S = diag([1/8;1;1;1;1;1;1;1;1;1;1;1;1]);
%S = eye(13);

%solve numerically via polynomial approximation
x0 = [7;0;0;0;0;0;0;0;0;1;0;0;0];
d = 13;

Dd = kron(D, eye(d));
u = repmat(S * x0, N+1, 1);

iSx = kron(eye(N+1), inv(S));

%solve by Newton iteration method
tic
for i =1:3
    xdot = f(u, u0, N, dynamics, S);
    Df = Jf(u, u0, N, dyn_jacobian, S);
    
    %compute step
    unew = (Dd * iSx - Df * iSx) \ -(xdot - (Dd * iSx) * u);
    u = u - unew;
end
toc

u = kron(eye(N+1), inv(S)) * u;
u_plot = reshape(u, 13, N+1);
plot(x,u_plot(1:3,:),'--');

%CHECK solutions
norm(Dd * u - f(u, u0, N, dynamics, S), inf)

norm(u_plot(10:13, 1))

end

function[xdot] = f(x,u, N, func, Scaling)
% u assumed to be constant
%scale variable
x = kron(eye(N+1), inv(Scaling)) * x;
x = reshape(x, 13, N+1);
%evaluate function at each point
xdot = [];
for i = 1:N
    eval = func({x(:,i), u});
    eval = full(eval{1});
    xdot = [xdot; eval];
end
xdot = [xdot;x(:,N+1)];
end


function[Df] = Jf(x, u, N, Jacobian, Scaling)
%evaluate system Jacobian
Df = zeros(13 * (N + 1));

d = 13; %state dimension
%scaling
x = kron(eye(N+1), inv(Scaling)) * x;
x = reshape(x, d, N+1);

%fill in Jacobian
for i = 0:N-1
    Jc = Jacobian({x(:,i+1), u});
    Df((d*i+1):(d*i+d),(d*i+1):(d*i+d)) = full(Jc{1});
end

end

