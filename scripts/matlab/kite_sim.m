function [NUM, FLOG, SYM] = kite_sim(aircraft, params)
%casadi based kite dynamics simulation
import casadi.*

%-------------------------
%Enviromental constants
%-------------------------
g = 9.80665; % gravitational acceleration [m/s2] [WGS84]
ro = 1.2985; % standart atmospheric density [kg/m3] [Standart Atmosphere 1976]

%---------------------------
%Glider geometric parameters
%---------------------------
b = aircraft.geometry.b;
c = aircraft.geometry.c;
AR = aircraft.geometry.AR;
S = aircraft.geometry.S;
lam = aircraft.geometry.lam;
St = aircraft.geometry.St; 
lt = aircraft.geometry.lt;
Sf = aircraft.geometry.Sf; 
lf = aircraft.geometry.lf;   
Xac = aircraft.geometry.Xac; 
Xcg = 0.031/c;               % Center of Gravity (CoG) wrt leading edge [1/c]
Vf = (Sf * lf) / (S * b);    % fin volume coefficient []
Vh = (lt * St) / (S * c);    % horizontal tail volume coefficient [] 

%---------------------------
%Mass and inertia parameters
%---------------------------
Mass = aircraft.inertia.mass;
Ixx = aircraft.inertia.Ixx;
Iyy = aircraft.inertia.Iyy; 
Izz = aircraft.inertia.Izz;
Ixz = aircraft.inertia.Ixz;

%-------------------------------
%Static aerodynamic coefficients
%-------------------------------
% All characteristics assumed linear
CL0 = aircraft.aerodynamic.CL0;
CL0_t = aircraft.aerodynamic.CL0_tail; 
CLa_tot = aircraft.aerodynamic.CLa_total;  
CLa_w = aircraft.aerodynamic.CLa_wing; 
CLa_t = aircraft.aerodynamic.CLa_tail;
e_o = aircraft.aerodynamic.e_oswald;
dw = CLa_tot / (pi * e_o * AR); % downwash acting at the tail []
CD0_tot = aircraft.aerodynamic.CD0_total; 
CD0_w = aircraft.aerodynamic.CD0_wing;
CD0_t = aircraft.aerodynamic.CD0_tail; 
CYb  = aircraft.aerodynamic.CYb;
CYb_vt = aircraft.aerodynamic.CYb_vtail;
Cm0 = aircraft.aerodynamic.Cm0;
Cma = aircraft.aerodynamic.Cma;
Cn0 = aircraft.aerodynamic.Cn0;
Cnb = aircraft.aerodynamic.Cnb; 
Cl0 = aircraft.aerodynamic.Cl0; 
Clb = aircraft.aerodynamic.Clb;

CLq = aircraft.aerodynamic.CLq;
Cmq = aircraft.aerodynamic.Cmq;
CYr = aircraft.aerodynamic.CYr; 
Cnr = aircraft.aerodynamic.Cnr; 
Clr = aircraft.aerodynamic.Clr;
CYp = aircraft.aerodynamic.CYp;
Clp = aircraft.aerodynamic.Clp; 
Cnp = aircraft.aerodynamic.Cnp; 

%------------------------------
%Aerodynamic effects of control
%------------------------------
CLde = aircraft.aerodynamic.CLde;
CYdr = aircraft.aerodynamic.CYdr;
Cmde = aircraft.aerodynamic.Cmde; 
Cndr = aircraft.aerodynamic.Cndr; 
Cldr = aircraft.aerodynamic.Cldr;
CDde = aircraft.aerodynamic.CDde; % (assume negligible)

CL_daoa = -2 * CLa_t * Vh * dw; % aoa-rate effect on lift (from Etkin) [] Stengel gives positive estimation !!!
Cm_daoa = -2 * CLa_t * Vh * (lt/c) * dw; %aoa-rate effect on pitch moment [] Stengel gives positive estimation !!!

%--------------------------
%State variables definition
%--------------------------
r = SX.sym('r', 3); % position of the CoG in the Inertial Reference Frame (IRF) [m]
q = SX.sym('q', 4); % body attitude relative to IRF [unit quaternion]
v = SX.sym('v', 3); % linear velocity of the CoG in the Body Reference Frame (BRF) [m/s]
w = SX.sym('w', 3); % glider angular velocity in BRF [rad/s]

%----------------------------
%Control variables definition
%----------------------------
%TODO: more detatailed propeller model will be added later
T = SX.sym('T'); % propeller propulsion : applies along X-axis in BRF [N]
dE = SX.sym('dE'); % elevator deflection [positive down] [rad]
dR = SX.sym('dR'); % rudder deflection [rad]
dA = SX.sym('dA'); % aileron deflection [reserved, but not used] [rad]
dF = SX.sym('dF'); % flaps deflection [reserved, but not used]

%Read algorithm and simulation parameters
%DEFAULTS
SIMULATION = 0;
USE_CVODES = 1;
VR_SIM = 0;
PLOT_RESULTS = 0;
X0 = [];
U0 = [];
T_START = 0.0;
T_FINISH = 1.0;
DT = 0.01;

if(isfield(params, 'simulation'))
    SIMULATION = params.simulation;
end

if(isfield(params, 'int_type'))
    if(isequal(params.int_type,'cvodes'))
        USE_CVODES = 1;
    else if (isequal(params.int_type, 'rk4'))
            USE_CVODES = 0;
        else
            fprintf(2, 'ERROR: Unknown integrator type: [%s] \n', params.int_type);
            error('Use "cvodes" or "rk4" instead');
        end
    end
end

%simulation time span
if(isfield(params, 't_span'))
    if(numel(params.t_span) ~= 3)
        error('ERROR: time span format should be: [start, finish, dt]');
    else
        T_START = params.t_span(1);
        T_FINISH = params.t_span(2);
        DT = params.t_span(3);
        if(T_START > T_FINISH)
            error('ERROR: T_START > T_FINISH');
        end
    end
end

if(SIMULATION)
    %VR visualisation
    if(isfield(params, 'vr'))
        VR_SIM = params.vr;
    end
    
    %initial condition and input
    if(isfield(params,'x0') && isfield(params,'u0'))
        X0 = params.x0;
        U0 = params.u0;
    else
        error('ERROR: initial conditions "x0" and control "u0" should be specified for simulation');
    end
    
    if(isfield(params,'plot'))
        PLOT_RESULTS = params.plot;
    end
end


%TODO: for now assume there is no wind, this functionality will be added later
V = L2(v);
V2 = v' * v;

ss = asin(v(2) / V); %side slip angle [rad] (v(3)/v(1)) // small angle assumption
aoa = atan2(v(3) , v(1));  %angle of attack definition [rad] (v(2)/L2(v))
dyn_press = 0.5 * ro * V2; %dynamic pressure

CD_w = CD0_w + (CL0 + CLa_w * aoa )^2 / (pi * e_o * AR); %wing drag coefficient
CD = CD0_tot + (CL0 + CLa_tot * aoa)^2 / (pi * e_o * AR); %total drag coefficient 

%-------------------------
%Dynamic Equations: Forces
%-------------------------
LIFT = (CL0 + CLa_tot * aoa) * dyn_press * S + ...
       (0.25 * CLq * c * S * ro) * V * w(2);
DRAG = CD * dyn_press * S;
SF = (CYb * ss + CYdr * dR) * dyn_press * S + ...
        0.25 * (CYr * w(3) + CYp * w(1)) * (b * ro * S) * V;

%Compute transformation betweeen WRF and BRF: qw_b
%qw_b = q(aoa) * q(-ss);
q_aoa = [cos(aoa/2); 0; sin(aoa/2); 0];
q_ss = [cos(-ss/2); 0; 0; sin(-ss/2)];

qw_b = quat_mul(q_aoa, q_ss);
qw_b_inv = quat_inv(qw_b);

%The following is a complete disaster
%Aerodynamic forces in BRF: Faer0_b = qw_b * [0; -DRAG; SF; -LIFT] * qw_b_inv
qF_tmp = quat_mul(qw_b_inv, [0; -DRAG; 0; -LIFT]);
qF_q = quat_mul(qF_tmp, qw_b);
Faero_b = qF_q(2:4);

Zde = (-CLde) * dE * dyn_press * S;
FdE_tmp = quat_mul(quat_inv(q_aoa), [0; 0; 0; Zde]);
FdE = quat_mul(FdE_tmp, q_aoa);
FdE = FdE(2:4);

Faero_b = Faero_b + FdE + [0; SF; 0];

%Gravity force in BRF
qG = quat_mul(quat_inv(q),[0;0;0;g]);
qG_q = quat_mul(qG, q);
G_b = qG_q(2:4);

%Propulsion force in BRF
T_b = [T;0;0];

%Tether force
%value: using smooth ramp approximation
Lt = 2.31; %tether length
d_ = L2(r);
%spring term
%Rv = ((d_ - Lt)) * Heaviside(d_ - Lt, 1);
Ks = 500 * Mass; %300 Ks: 800;
Kd = 10 * Mass; %30 Kd: 500
Rv = ((d_ - Lt));
Rs = -Rv * (r / d_);
%damping term
qvi = quat_mul(q, [0; v]);
qvi_q = quat_mul(qvi, quat_inv(q));
vi = qvi_q(2:4);
Rd = (-r / d_) * (r' * vi) / d_;
R = ( Ks * Rs + Kd * Rd) * Heaviside(d_ - Lt, 1);

%BRF
qR = quat_mul(quat_inv(q),[0;R]);
qR_q = quat_mul(qR, q);
R_b = qR_q(2:4);


%Total external forces devided by glider's mass (linear acceleration)
v_dot = (Faero_b + T_b  + R_b)/Mass + G_b  - cross(w,v);

%-------------------------
%Dynamic Equation: Moments
%-------------------------
%Rolling Aerodynamic Moment
L = (Cl0 + Clb * ss + Cldr * dR) * dyn_press * S * b + ...
     (Clr * w(3) + Clp * w(1)) * (0.25 * ro * b^2 * S) * V;


%Pitching Aerodynamic Moment
M = (Cm0 + Cma * aoa  + Cmde * dE) * dyn_press * S * c + ...
     Cmq * (0.25 * S * c^2 * ro) * w(2) * V;

%Yawing Aerodynamic Moment
N = (Cn0 + Cnb * ss + Cndr * dR) * dyn_press * S * b + ...
     (Cnp * w(1) + Cnr * w(3)) * (0.25 * S * b^2 * ro) * V;
 
%Aircraft Inertia Matrix 
J = [Ixx 0 Ixz; 0 Iyy 00; Ixz 0 Izz];

%Angular motion equationin BRF
Maero = [L; M; N];
%Moments rotation SRF -> BRF
T_tmp = quat_mul(quat_inv(q_aoa), [0; Maero]);
Trot = quat_mul(T_tmp, q_aoa);
Maero = Trot(2:4);

w_dot = inv(J) * (Maero - cross(w, J * w));

%-----------------------------
%Kinematic Equations: Position
%-----------------------------
% Aircraft position in the IRF
qv = quat_mul(q, [0; v]);
qv_q = quat_mul(qv, quat_inv(q));
r_dot = qv_q(2:4);

%-----------------------------
%Kinematic Equations: Attitude
%-----------------------------
%Aircraft attitude wrt to IRF 
q_dot = 0.5 * quat_mul(q, [0;w]);

%Now we completely define dynamics of the aircraft
state = [v; w; r; q];
control = [T; dE; dR];
dynamics = [v_dot; w_dot; r_dot; q_dot];

dyn_func = Function('dynamics', {state, control}, {dynamics});

%compute dynamics state Jacobian
d_jacobian = dynamics.jacobian(state);
dyn_jac = Function('dyn_jacobian',{state, control}, {d_jacobian});

%define RK4 integrator scheme
X = SX.sym('X',13);
U = SX.sym('U',3);
dT = SX.sym('dT');

%get symbolic expression for an integrator
integrator = RK4_sym(X, U, dyn_func, dT);
RK4_INT = Function('RK4', {X,U,dT},{integrator});

%get Jacobian of the RK4 Integrator
integrator_jacobian = integrator.jacobian(X);
rk4_jacobian = Function('RK4_JACOBIAN', {X, U, dT}, {integrator_jacobian});

h = DT;
x0 = X0;
u0 = U0;

%CVODES integrator
ode = struct('x',state, 'p',control, 'ode',dynamics);
opts = struct('tf', h);
CVODES_INT = casadi.integrator('CVODES_INT','cvodes', ode, opts);

aoa_func = Function('AoA',{v},{aoa});
ss_func = Function('SS',{v},{ss});

%return integrator function
if(USE_CVODES)
    NUM.INT = CVODES_INT;
else
    NUM.INT = RK4_INT;
end

%Symbolic expression for the model
SYM.INT = integrator;
SYM.DYNAMICS = dynamics;
SYM.DYN_JAC = d_jacobian;

SYM.state = state;
SYM.control = control;

SYM.X = X;
SYM.U = U;
SYM.dT = dT;

NUM.DYN_FUNC = dyn_func;
NUM.DYN_JACOBIAN = dyn_jac;
NUM.RK4_JACOBIAN = rk4_jacobian;

if(VR_SIM)
    %create simulation scene
    [scene] = CreateScene('BackView',1);
end

%-----------------------------
%Open-Loop Flight Simulation
%-----------------------------

%simulation & visualisation
FLOG = [];
if (SIMULATION)
    t = T_START;
    log = [];
    flight_angles = [];
    res = [];
    while t <= T_FINISH
        %simulation loop
        %logging data
        log = [log; x0', t];
        
        aoa_t = aoa_func({x0(1:3)});
        full(aoa_t{1});
        ss_t = ss_func({x0(1:3)});
        full(ss_t{1});
        flight_angles = [flight_angles; full(aoa_t{1}) full(ss_t{1})];
        
        if (USE_CVODES)
            out = dyn_func({x0,u0});
            dyn =full(out{1});
            
            out = CVODES_INT(struct('x0',x0, 'p',u0));
            res = full(out.xf);
        else
            out = RK4_INT({x0, u0, h});
            res = full(out{1});
        end
        
        t = t + h;
        x0 = res;
        
        if(VR_SIM)
            UpdateScene(scene, res(7:9), res(10:13),t);
        end
        
    end
    
    %return flight telemetry
    FLOG = log;
    
    if(PLOT_RESULTS)
        %trajectory 
        figure
        plot3(log(:,7), log(:,8), -log(:,9));
        axis equal
        grid on
        
        %velocities in BRF
        figure('Name','Linear Velocities in BRF','NumberTitle','off')
        plot(log(:,end), log(:,1:3))
        grid on
        
        %angle of attack
        figure('Name','Angle of Attack & Sideslip Angle','NumberTitle','off')
        %plot(rad2deg(flight_angles));
        plot(log(:,end), flight_angles);
        grid on
        
        %attitude
        figure('Name','Attitude (EulerAngles)','NumberTitle','off')
        plot(log(:,end), quat2eul(log(:, 10:13),'ZYX'))
        grid on
        
        %angular rates
        figure('Name','Angular rates in BRF','NumberTitle','off')
        plot(log(:,end), log(:,4:6))
        grid on
    end
    
end

end

%Casadi helper functions

%L2-norm of a vector
function [L2] = L2(x)
if(size(x,1) >= size(x,2))
    L2 = sqrt(x' * x);
else
    L2 = sqrt(x * x');
end

end

%quaternion multiplication
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

%Runge-Kutta 4th-order symbolic integration method
function [x_h] = RK4_sym(x, u, f, h)
% this integrator is problem specific
k1 = f({x,u});
k2 = f({x + 0.5 * h * k1{1},u});
k3 = f({x + 0.5 * h * k2{1},u});
k4 = f({x + h * k3{1},u});

x_h = x + (h/6) * (k1{1} + 2*k2{1} + 2*k3{1} + k4{1});
end

%Heaviside smooth approximation
function[H] = Heaviside(x,k)
H = k / (1 + exp(-4*x));
end

function [scene] = CreateScene(viewpoint, is_recording)
%create 3D scene
scene = struct;

scene.world = vrworld('world.wrl','new');
open(scene.world);
fig = vrfigure(scene.world);
% go to a viewpoint suitable for user navigation
set(fig, 'Viewpoint', viewpoint);

% get the manipulated airplane node
scene.position = vrnode(scene.world, 'GliderTranslational');
scene.attitude = vrnode(scene.world, 'GliderRotational');

if(is_recording)
    set(scene.world, 'RecordMode', 'manual');
    set(fig, 'Record2DFileName', 'flight_sim.avi');
    set(fig, 'Record2D', 'on');
    set(fig, 'Record2DCompressQuality', 100);
end
end

function [] = UpdateScene(scene, position, attitude, time)
%frame transformations with VR simulator
q_sim = eul2quat([ 0 0 pi/2]);

Xi = position;
Qi = attitude;

%apply simulation and draw
x_sim = quatmul(quatmul(q_sim, [0 Xi']),quatinv(q_sim));
q_disp = quatmul(quatmul(q_sim, Qi'),quatinv(q_sim));

scene.position.translation = x_sim(2:4);
scene.attitude.rotation = quat2axang(q_disp);

set(scene.world, 'Time', time);

vrdrawnow;
end