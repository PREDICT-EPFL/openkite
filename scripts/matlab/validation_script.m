%validation script
clc
%raw_data_filename = '../flight1.log';
raw_data_filename = '../flight_banked.log';
%raw_data_filename = '../flight_elev_up.log';
%raw_data_filename = '../flight_rudder_neg.log';
%raw_data_filename = '../flight_thrust.log';

%load the raw data
data = load(raw_data_filename);

%substract interesting data
idx_start = 1300; %thrust 2140; %rudder neg 880;  %elev_up 1200; banked 1300;  %level flight 780;
idx_finish = 1550; %thrust 2365; %rudder neg 1210; %elev_up 1820; banked 1550; %level flight 1150;

idata = data(idx_start:idx_finish, :);

%transform data to the IRF 
transformation = eul2quat([0 0 -pi]);
brf_offset = [-0.09; 0; -0.02];

irf_data = optitrack2world(idata, transformation, brf_offset);

%cook the telemetry data
exp_telemetry = [];
identity = [1 0 0 0];
vw = [];
dt = 0;
dt_log = [];

rdot_log = [];


for i = 1:size(irf_data,1)-1
    %compute state derivatives using finite differences
    dt = irf_data(i+1,1) - irf_data(i,1);
    dt_log = [dt_log; dt];
    rdot = (irf_data(i+1, 2:4) - irf_data(i,2:4)) ./ dt;
    att = irf_data(i,5:8);
    att_inv = quatinv(irf_data(i,5:8));
    v_body = quatmul(quatmul(att_inv, [0 rdot]), quatinv(att_inv));
    %linear velocity in the BRF
    v_body = v_body(2:4);
    
    %angular rates in BRF
    att2 = irf_data(i+1,5:8);
    dq = quatmul(att_inv, att2);
    
    qdot = dq - identity;
    q_w = 2 * (dq - identity) ./ dt;
    w_body = q_w(2:4);
    
    vw = [vw; v_body w_body];
end

%smoothen data with median filter
m_filter_order = 1;
smooth_vw = medfilt1(vw(:,4:6), m_filter_order);
dt = mean(dt_log);

%compile telemetry log
exp_telemetry = [exp_telemetry vw(:,1:3) smooth_vw];
%append position, attitude and time
exp_telemetry = [exp_telemetry irf_data(1:end-1, 2:8)];
exp_telemetry = [exp_telemetry irf_data(1:end-1, 1)];

%playback the experimental flight
%flight_playback(exp_telemetry)

%setup simulation
aircraft = ReadYaml('umx_radian.yaml');
parameters.simulation = 1;
parameters.plot = 0;
parameters.vr = 0;
parameters.int_type = 'rk4';

start_sample = 30; %15 elev_up; 30 banked; %50  level;
s_time = exp_telemetry(start_sample,end);
f_time = exp_telemetry(end,end);
parameters.t_span = [s_time, f_time, dt];
parameters.u0 = zeros(3,1);
parameters.x0 = exp_telemetry(start_sample,1:13)';
%control
parameters.u0 = [0;0;0];

disp(aircraft)
disp(parameters)

%run simulation
%[int, flog] = flight_sim(aircraft, parameters);
[num, flog, sym] = kite_sim(aircraft, parameters);

%check state reconstruction with forward integration
FORWARD_INT = 0;
if(FORWARD_INT)
    x0_FI = exp_telemetry(start_sample,7:13);
    log_FI = [x0_FI];
    RB_INT = rigid_body(dt);
    for i = start_sample:length(exp_telemetry)-1
        vb = exp_telemetry(i,1:3);
        wb = exp_telemetry(i,4:6);
        out = RB_INT(struct('x0', x0_FI, 'p', [vb';wb']));
        x0_FI = full(out.xf);
        log_FI = [log_FI; x0_FI'];
    end
    
    size(flog)
    size(log_FI)
end

%Extended Kalman Filter test
V = diag([0.005;0.005;0.005; 0.001;0.001;0.001;0.001].^2);
%W = zeros(13);                                        % model noise
x_est = exp_telemetry(start_sample,1:13)';            % state estimation 
e_v = [0.001; 0.001; 0.001];
e_w = [0.5; 0.5; 0.05];
e_r = [0.001; 0.001; 0.001];
e_q = [0.001; 0.001; 0.001; 0.001];
P_est = 100 * diag([e_v; e_w; e_r; e_q].^2);                 % estimation covariance
W = diag([e_v; e_w; e_r; e_q].^2);

model.num = num;
model.sym = sym;

estim = [x_est'];

dyn_func = num.DYN_FUNC;

%for k = start_sample:(size(exp_telemetry,1)-1)
for k = start_sample:100
    dt = exp_telemetry(k+1,end) - exp_telemetry(k,end);
    %get estimation
    u = parameters.u0;
    z_m = exp_telemetry(k+1,7:13);

    [x_est, P_est] = kiteEKF(x_est, P_est, u, model, dt, W, V, z_m');
    estim = [estim; x_est'];
end

%flight_playback(flog);
%check state reconstruction with forward integration
FORWARD_INT = 1;
if(FORWARD_INT)
    x0_FI = estim(1,7:13);
    log_FI = [x0_FI];
    %dt = mean(dt_log);
    RB_INT = rigid_body(dt);
    for i = 1:length(estim)-1
        vb = estim(i,1:3);
        wb = estim(i,4:6);
        out = RB_INT(struct('x0', x0_FI, 'p', [vb';wb']));
        x0_FI = full(out.xf);
        log_FI = [log_FI; x0_FI'];
    end
    
    size(flog)
    size(log_FI)
end

norms = sqrt(sum(exp_telemetry(start_sample:end,1:3).^2,2));
aoa = atan2(exp_telemetry(start_sample:end,3), exp_telemetry(start_sample:end,1));
ss = asin(exp_telemetry(start_sample:end,2)./norms);
%plot(exp_telemetry(start_sample:end,end), aoa);
%hold on
%plot(exp_telemetry(start_sample:end,end), ss);
%grid on

%plot results
figure('Name','Linear Velocities','NumberTitle','off')
plot(exp_telemetry(:,end), exp_telemetry(:,1:3),'--')
grid on
hold on
plot(flog(:,end), flog(:,1),'b', flog(:,end), flog(:,2),'r', flog(:,end), flog(:,3),'g')
plot(flog(1:size(estim,1),end), estim(:, 1:3),'o');

figure('Name','Angular Velocities','NumberTitle','off')
plot(exp_telemetry(:,end), exp_telemetry(:,4:6),'--')
grid on
hold on
plot(flog(:,end), flog(:,4),'b', flog(:,end), flog(:,5),'r', flog(:,end), flog(:,6), 'g')
plot(flog(1:size(estim,1),end), estim(:, 4:6),'o');

%trajectory
figure('Name','Trajectory')
plot3(exp_telemetry(:,7), exp_telemetry(:,8), -exp_telemetry(:,9));
grid on
axis equal
hold on
plot3(flog(:,7), flog(:,8), -flog(:,9))
plot3(log_FI(:,1), log_FI(:,2), -log_FI(:,3))
plot3(estim(:,7), estim(:,8), -estim(:,9))


%attitude
rpy_exp = quat2eul(exp_telemetry(:,10:13));
rpy_sim = quat2eul(flog(:,10:13));
rpy_est = quat2eul(log_FI(:,4:7));

figure('Name','Attitude (Quaternion)','NumberTitle','off')
plot(exp_telemetry(:,end), rpy_exp,'--')
grid on
hold on
plot(flog(:,end), rpy_sim(:,1),'b',flog(:,end), rpy_sim(:,2), 'r', flog(:,end), rpy_sim(:,3), 'g');
plot(flog(1:size(rpy_est,1),end), rpy_est,'o');

