function[f_log, c_log] = kite_collocation()
%Chebyshev collocation based path following for a glider
import casadi.*

aircraft = ReadYaml('umx_radian.yaml');
parameters.simulation = 0;
parameters.plot = 0;
parameters.vr = 0;
parameters.int_type = 'cvodes';

disp(aircraft)
disp(parameters)

%get aircraft model
[num, log, sym] = kite_sim(aircraft, parameters);

%get dynamics function and state Jacobian
dynamics = sym.DYNAMICS;
X = sym.state;
U = sym.control;

%state and control dimensionality
n = 15;
m = 4;

%Order of polynomial interpolation
N = 12;

%define augmented sdynamics of path parameter
V = SX.sym('V', 2);
Uv = SX.sym('Uv');
Av = [0 1; 0 0];
Bv = [0;1];

%parameter dynamics
p_dynamics = Av * V + Bv * Uv;

%augmented system
aug_state = [X; V];
aug_control = [U; Uv];
aug_dynamics = [dynamics; p_dynamics];

%evaluate augmented dynamics
aug_dynamo = Function('AUG_DYNAMO', {aug_state, aug_control}, {aug_dynamics});

%NLP formulation
%set up control constraints --------------------------------------------- %

Tmin = 0;                 %[N]
Tmax = 0.3;               %[N]
elev_min = deg2rad(-8);   %[rad]
elev_max = deg2rad(8);    %[rad]
rudder_min = deg2rad(-8); %[rad]
rudder_max = deg2rad(8);  %[rad]
uz_min = -10.0;           %[m/s2]
uz_max = 10.0;            %[m/s2]

%state constraints ------------------------------------------------------ %
wx_min = -2*pi;  %[rad/s]
wx_max = 2*pi;   %[rad/s]
wy_min = -2*pi;  %[rad/s]
wy_max = 2*pi;   %[rad/s]
wz_min = -2*pi;  %[rad/s]
wz_max = 2*pi;   %[rad/s]

vx_max = 12.0;  %[m/s]
vx_min = 0.5;   %[m/s]
vy_max = 5.0;   %[m/s]
vy_min = -5.0;   %[m/s]
       
%unscaled state constraints
z_min = [vx_min; vy_min; -inf; wx_min; wy_min; wz_min; ...
         -inf; -inf;-inf; -1; -1; -1; -1; -inf; -inf];
      
z_max = [vx_max;  vy_max; inf; wx_max; wy_max; wz_max; ... 
           inf; inf; inf; 1; 1; 1; 1; inf; inf];
              
%unscaled input constraints
u_min = [Tmin; elev_min; rudder_min; uz_min];
u_max = [Tmax; elev_max; rudder_max; uz_max];

%weight matrices for the cost function ---------------------------------- %
Q = diag([1e4; 1e4; 5e3]);      %error penalty
R = 1e-4 * diag([1; 1; 1; 1]);  %control derivative penalty
W = 1e-3;                       %velocity profile following 
Wq = diag([1000,1000]);

%reference velocity for the parameter
vel_ref = 0.05 ;    %[1/s] ? change to kite airspeed ?
rsphere = 5;  %[m] max elongation of a tether

%scaling matrices --------------------------------------------------------%
%state
Sx = diag([1/vx_max; 1/vy_min; 1; 1/wx_max; 1/wy_max; 1/wz_max; ... 
           1/rsphere;1/rsphere;1/rsphere; 1; 1; 1; 1; 1/(2*pi); 1/vx_max]);
invSx = inv(Sx);
%control
Su = diag(abs([1/Tmax; 1/elev_max; 1/rudder_max; 1/uz_max]));
invSu = inv(Su);

%scaled constraints ------------------------------------------------------%
z_max_scaled = Sx * z_max;
z_min_scaled = Sx * z_min;
u_max_scaled = Su * u_max;
u_min_scaled = Su * u_min;

%define optimisation variables
z = SX.sym('z', n, N+1);
z_u = SX.sym('z_u', m, N);

%allocate vectors to store constraints
vecx = {};
g = {};

%state and control constraints
lbx = [];
ubx = [];

%nonlinear state constraints
lbg = [];
ubg = [];

%objective function
obj = 0;

%Geometric path definition -----------------------------------------------%
% we will fly a simple horizontal circle
radius = aircraft.tether_length + 0.1;   %[m]
altitude = 0.0;                          %[m]
theta = SX.sym('theta');
path_xyz = [radius * cos(theta); radius * sin(theta); altitude];
path_fun = Function('path', {theta}, {path_xyz});

%scale
path_fun_scaled = Function('path', {theta}, {Sx(7:9,7:9) * path_xyz});

%---- NLP formulation --------------------%
%---- Lagrangian term
pos = z(7:9, 2:end-1);
pos = pos(:);

%scaled
path = f_vec(Sx(14,14) * z(14,2:end-1), path_fun_scaled);

%path following error
err = pos - path;
Qn = kron(eye(N-1), Q);

%reference velocity following
err_v = z(15, 2:end-1)' - repmat(vel_ref, N-1, 1);
Qw = W * eye(N-1, N-1);

L = err' * Qn * err + err_v' * Qw * err_v;

%Regularization term
%Ru = kron(eye(N), R);
%L = L + z_u(:)' * Ru * z_u(:);

%---- Mayer term ------------------------%
%path = path_fun({z(14,1)});
%scaled
path = path_fun_scaled({Sx(14,14) * z(14,1)});
err_n = z(7:9,1) - path{1};
err_n = [err_n; z(15, 1) - vel_ref];

T = err_n' * (10 * blkdiag(Q, W)) * err_n;

%objective function
obj = L + T;

%----- set equality constraints ---------% 
%differentiation matrix
[x, D] = cheb(N);
D(end,:) = zeros(N+1,1);
D(end, end) = 1;
Dn = kron(D, eye(n));

%scaling
iSx = kron(eye(N+1), invSx);
iSu = kron(eye(N), invSu);

%initial condition
z0 = z(:,end);

F = f(invSx * z, invSu * z_u, z0, N, aug_dynamo);
G = Dn * iSx * z(:) - F;
% G = Dn * z(:) - kron(eye(N+1), Sx) * F; ? try that instead

g = {g{:}, G(1:(n * N))};
lbg = [lbg, zeros((N) * n,1)'];
ubg = [ubg, zeros((N) * n,1)'];

%----- set inequality (box) constraints -----------------%
%state
lbx = [lbx; repmat(z_min_scaled, N+1, 1)];
ubx = [ubx; repmat(z_max_scaled, N+1, 1)];

%control
lbx = [lbx; repmat(u_min, N, 1)];
ubx = [ubx; repmat(u_max, N, 1)];

%----- set decision variable ---------------------------%
vecx = {vecx{:}, z(:), z_u(:)};

%formulate NLP
nlp = struct('x',vertcat(vecx{:}), 'f', obj, 'g', vertcat(g{:}));

%solver = nlpsol('solver','sqpmethod',nlp, struct('qpsol', 'qpoases'));
options = struct('ipopt', struct('linear_solver','ma97'));
solver = nlpsol('solver', 'ipopt', nlp, options);

%initial_guess
feasible_initial_state = Sx * [1.5; 0; 0; 0; 0; 0;
                               0; 1.0; 0; 1.0; 0; 0; 0; pi/2; 0];
%scaled
feasible_initial_control = Su * [0.1;0;0;0];

args = struct;
args.lbx = lbx;
args.ubx = ubx;
args.lbg = lbg;
args.ubg = ubg;
args.x0 = repmat(feasible_initial_state, N+1, 1);
args.x0 = [args.x0; repmat(feasible_initial_control, N, 1)];
%update state constraints
args.lbx((n * N)+1:n*(N+1)) = feasible_initial_state;
args.ubx((n * N)+1:n*(N+1)) = feasible_initial_state;

res = solver(args);
x = full( res.x );

opt_x = x(1:(n * (N+1)));
opt_x = reshape(opt_x, n, N+1);

%get control
opt_u = x(((n * (N+1)+1):end));
opt_u = reshape(opt_u, m, N);

% ---------------- Closed Loop Simulation --------------------- %

sim_time = 2.0;          %[s]
time = 0;
h = 0.02;
time_line = [];
f_log     = [];
ctl_log   = [];
delay_log = [];

% states-control subplots code
figure
h1 = subplot(6,3,1);
h2 = subplot(6,3,2);
h3 = subplot(6,3,3);
h4 = subplot(6,3,4);
h5 = subplot(6,3,5);
h6 = subplot(6,3,6);
h7 = subplot(6,3,7);
h8 = subplot(6,3,8);
h9 = subplot(6,3,9);
h10 = subplot(6,3,10);
h11 = subplot(6,3,11);
h12 = subplot(6,3,12);
h13 = subplot(6,3,13);
h14 = subplot(6,3,14);
h15 = subplot(6,3,15);
h16 = subplot(6,3,16);
h17 = subplot(6,3,17);
h18 = subplot(6,3,18);
grid on

% trajectory visualization
traj = figure;
figure(traj);
path_points = compute_path([0:0.1:2*pi], path_fun);
plot3(path_points(1,:), path_points(2,:), path_points(3,:),'r');
grid on
hold on
axis equal

%CVODES integrator for augmented system
ode = struct('x',aug_state, 'p',aug_control, 'ode',aug_dynamics);
opts = struct('tf', h);
CVODES_INT = casadi.integrator('CVODES_INT','cvodes', ode, opts);

while time <= sim_time
    %integrate flight dynamics
    %rescale and integrate
    x_ = invSx * opt_x(:,end);
    u_ = invSu * opt_u(:,end);
    x_next = CVODES_INT(struct('x0', x_, 'p', u_));
    x_next = full(x_next.xf);
    x_next(10:13) = x_next(10:13) / norm(x_next(10:13));
    
    %update time step
    time = time + h;
    
    %put telemetry in the log
    time_line = [time_line; time];
    f_log = [f_log; x_next'];
    
    %scale back to optimization problem
    x_next = Sx * x_next;
    
    %update constraints and resolve NLP
    args.lbx((n * N)+1:n*(N+1)) = x_next;
    args.ubx((n * N)+1:n*(N+1)) = x_next;
    
    %initial guess
    args.x0 = res.x;
    
    %solve
    tic
    res = solver(args);
    delay = toc;
    delay_log = [delay_log; delay];
    
    %separate variables
    x = full( res.x );
    opt_x = x(1:(n * (N+1)));
    opt_x = reshape(opt_x, n, N+1);
    opt_u = x(((n * (N+1)+1):end));
    opt_u = reshape(opt_u, m, N);
    
    %update control log
    %scaled
    u_log = invSu * opt_u(:,end);
    ctl_log = [ctl_log; u_log'];
    
    eul = quat2eul(f_log(:,10:13));
    
    plot(h1, time_line, f_log(:,1)); grid(h1,'on');
    plot(h2, time_line, f_log(:,2)); grid(h2,'on');
    plot(h3, time_line, f_log(:,3)); grid(h3,'on');
    plot(h4, time_line, f_log(:,4),'b',time_line, repmat(wx_max,length(time_line),1),'r', time_line, repmat(wx_min,length(time_line),1),'r'); grid(h4,'on');
    plot(h5, time_line, f_log(:,5),'b',time_line, repmat(wy_max,length(time_line),1),'r', time_line, repmat(wy_min,length(time_line),1),'r'); grid(h5,'on');
    plot(h6, time_line, f_log(:,6),'b',time_line, repmat(wz_max,length(time_line),1),'r', time_line, repmat(wz_min,length(time_line),1),'r'); grid(h6,'on');
    plot(h7, time_line, f_log(:,7)); grid(h7,'on');
    plot(h8, time_line, f_log(:,8)); grid(h8,'on');
    plot(h9, time_line, -f_log(:,9),'b', time_line, repmat(0.5,length(time_line),1),'r'); grid(h9,'on');
    plot(h10, time_line, eul(:,3)); grid(h10,'on');
    plot(h11, time_line, eul(:,2)); grid(h11,'on');
    plot(h12, time_line, eul(:,1)); grid(h12,'on');
    stairs(h13, time_line, ctl_log(:,1)); grid(h13,'on');
    stairs(h14, time_line, ctl_log(:,2)); grid(h14,'on');
    stairs(h15, time_line, ctl_log(:,3)); grid(h15,'on');
    plot(h16, time_line, f_log(:,14)); grid(h16,'on');
    plot(h17, time_line, f_log(:,15)); grid(h17,'on');
    plot(h18, time_line, ctl_log(:,4)); grid(h18,'on'); 
    
    %draw results
    plot_x = world2matlab(invSx(7:9,7:9) * opt_x(7:9,:));
    plot3(plot_x(1,:), plot_x(2,:), plot_x(3,:),'b');
    plot3(plot_x(1,end), plot_x(2,end), plot_x(3,end),'b*');
    path_points = compute_path(Sx(14,14) * opt_x(14,:), path_fun);
    plot3(path_points(1,:), path_points(2,:), path_points(3,:),'mo');
    
    drawnow;
end

figure
plot(delay_log, '*');
grid on;

end


function[xdot] = f(x, u, x0, N, func)
%evaluate function at each point
xdot = [];
for i = 1:N
    eval = func({x(:,i), u(:,i)});
    xdot = [xdot; eval{1}];
end
xdot = [xdot; x0];
end

function[out] = f_vec(x, fun)
out = [];
for i = 1:length(x)
    res = fun({x(i)});
    out = [out; res{1}];
end
end

%transform kite trajectory to matlab plot reference frame
function[rvis] = world2matlab(trajectory)
traj_size = size(trajectory,2);
rvis = zeros(3, traj_size);
transform = eul2quat([0, deg2rad(180), 0]);
for i = 1:traj_size
    pos = trajectory(:,i);
    res = quatmul(quatmul(transform, [0 pos']), quatinv(transform)); 
    rvis(:,i) = res(2:4);
end
end

%compute path points
function[path] = compute_path(param, path_fun)
path = zeros(3, length(param));
for i = 1 : length(param)
    res = path_fun({param(i)});
    path(:,i) = full(res{1});
end
end
