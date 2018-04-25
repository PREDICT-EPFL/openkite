%kite validation script
clc
control_filename = 'kite_control.log';
pose_filename = 'kite_pose.log';
%state_filename = 'kite_state.log';

%load the raw data
control = load(control_filename);
pose = load(pose_filename);
%state = load(state_filename);

%substract pivoting point
%pivot = [ -45.7; -70.3; 2824.2] / 1000;
%pivot = [0.189731; -0.04971; 2.8603411];
 pivot = [-0.1564; -0.2368; 2.8547];
for i = 1:length(pose)
    pose(i, 2:4) = pose(i, 2:4) - pivot';
end

plot3(pose(:,2), pose(:,3), pose(:,4))
grid on
axis equal

figure
dist = sqrt(sum(abs(pose(:,2:4)).^2,2));
plot(dist)
grid on

%%
%fit the sphere
%center = pivot;
%radius = 2.45;
%[center, radius, error] = sphere_fit_3d(pose(4000:10000, 2:4), center, radius)

%%
%transform pwm to control signal
max_throttle = 0.3;
max_rudder = deg2rad(8);
max_elevator = deg2rad(8);
throttle = (control(:, 2) - min(control(:,2))) * max_throttle / ...
           (max(control(:,2)) - min(control(:, 2)));

mid = 1470;       
rudder = (control(:,3) - mid) * max_rudder / ...
          (max(control(:, 3)) - min(control(:, 3)));

elevator = (control(:,4) - mid) * max_rudder / ...
          (max(control(:, 4)) - min(control(:, 4)));
      
control = [control(:, 1), throttle, elevator, rudder];

%transform data to the IRF 
transformation = eul2quat([0 0 -pi]);
brf_offset = [-0.09; 0; -0.02];
pose_irf = optitrack2world(pose, transformation, brf_offset);

%% translate and substract pivot point
% state check
pivot_irf = quatmul(quatmul(transformation, [0 pivot']), quatinv(transformation));
pivot_irf = pivot_irf(2:4);

for i = 1:length(state)
    state(i, 8:10) = state(i, 8:10) - pivot_irf;
end

%%
plot3(pose_irf(:,2), pose_irf(:,3), -pose_irf(:,4))
grid on
hold on
axis equal
plot3(state(:,8), state(:,9), -state(:,10))

%%
%find common start point
start = min([pose_irf(1,1), control(1,1)]);
control(:,1) = control(:,1) - start;
pose_irf(:,1) = pose_irf(:,1) - start;

%cook the telemetry data
exp_telemetry = [];
identity = [1 0 0 0];
vw = [];
dt_log = [];
rdot_log = [];

for i = 1:size(pose_irf, 1)-1
    %compute state derivatives using finite differences
    dt = pose_irf(i+1,1) - pose_irf(i,1);
    dt_log = [dt_log; dt];
    rdot = (pose_irf(i+1, 2:4) - pose_irf(i,2:4)) ./ dt;
    att = pose_irf(i, 5:8);
    att_inv = quatinv(pose_irf(i, 5:8));
    v_body = quatmul(quatmul(att_inv, [0 rdot]), quatinv(att_inv));
    %linear velocity in the BRF
    v_body = v_body(2:4);
    
    %angular rates in BRF
    att2 = pose_irf(i+1,5:8);
    dq = quatmul(att_inv, att2);
    
    qdot = dq - identity;
    q_w = 2 * (dq - identity) ./ dt;
    w_body = q_w(2:4);
    
    vw = [vw; v_body w_body];
end

%smoothen data with median filter
m_filter_order = 3;
smooth_w = medfilt1(vw(:,4:6), m_filter_order);
smooth_v = medfilt1(vw(:,1:3), m_filter_order);
dt = mean(dt_log);

%compile telemetry log
exp_telemetry = [exp_telemetry smooth_v smooth_w];
%append position, attitude and time
exp_telemetry = [exp_telemetry pose_irf(1:end-1, 2:8)];
exp_telemetry = [exp_telemetry pose_irf(1:end-1, 1)];

plot(exp_telemetry(:, 10:13))
grid on
%flight_playback(exp_telemetry)

%%
%create kite model
aircraft = ReadYaml('umx_radian.yaml');
parameters.simulation = 0;
parameters.plot = 0;
parameters.vr = 0;
parameters.int_type = 'cvodes';

start_sample = 10;
s_time = exp_telemetry(start_sample,end);
f_time = exp_telemetry(end,end);
parameters.t_span = [s_time, f_time, dt];
parameters.x0 = exp_telemetry(start_sample,1:13)';
parameters.u0 = [0;0;0];

[num, flog, sym] = kite_sim(aircraft, parameters);

%%

%check state reconstruction with forward integration
FORWARD_INT = 1;
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
    
    size(log_FI)
    save('forward_simulation_rb.mat', 'log_FI')
else
    load('forward_simulation_rb.mat');
end

%% 
% print forward integration results
plot3(exp_telemetry(:,7), exp_telemetry(:,8),-exp_telemetry(:,9))
hold on
plot3(log_FI(:,1), log_FI(:,2), -log_FI(:,3))
axis equal

%%

%Extended Kalman Filter test
V = diag([0.01;0.01;0.01; 0.0001;0.005;0.005;0.005].^2);
%W = zeros(13);                                        % model noise
x_est = exp_telemetry(start_sample,1:13)';            % state estimation
e_v = [0.5; 0.5; 0.5];             %[0.05; 0.5; 0.1];
e_w = [0.5; 0.5; 0.5];             %[0.1 0.1 0.5]
e_r = [0.5; 0.1; 0.1];             %[0.25; 0.25; 0.05];
e_q = [0.1; 0.05; 0.05; 0.05];    %[0.3; 0.05; 0.05; 0.06];
P_est = 10 * diag([e_v; e_w; e_r; e_q].^2);                 % estimation covariance
W = diag([e_v; e_w; e_r; e_q].^2);

model.num = num;
model.sym = sym;
model.num.DRE = dre(13, dt);

estim = [x_est'];

dyn_func = num.DYN_FUNC;
c_time = control(:,1);
points = 3000;

p_norms = [];
W_est = W;
west_log = diag(W_est);
%W = W_est;
res_log = [];
ctl_log = [];
delays = [];
x_sim = x_est;

rb_model = rigid_body(dt);
rb_model.num.DRE = model.num.DRE;

for k = start_sample:start_sample + points - 1
    dt = exp_telemetry(k+1,end) - exp_telemetry(k,end);
    delays = [delays;dt];
    %get corresponding control
    ctl_idx = length(c_time(c_time <= exp_telemetry(k, end)));
    if (isequal(ctl_idx, 0))
        u = [0,0,0];
    else
        u = control(ctl_idx, 2:4);
    end
    ctl_log = [ctl_log, u'];
    %get observation
    z_m = exp_telemetry(k+1,7:13);

    %propagate
    %x_sim(2) = x_est(2);
    %out = model.num.INT(struct('x0',x_est, 'p',u));
    %x_est = full(out.xf);
    [x_est, P_est, W_est, res] = kiteEKF(x_est, P_est, u, rb_model, dt, W, V, z_m', parameters.int_type);
    estim = [estim; x_est'];
    %west_log = [west_log, diag(W_est)];
    res_log = [res_log, res];
    p_norms = [p_norms, norm(real(P_est), 2)];
end

%%
idx_plot = start_sample:(start_sample + points);
figure
plot3(exp_telemetry(idx_plot, 7), exp_telemetry(idx_plot, 8), -exp_telemetry(idx_plot, 9))
%plot(exp_telemetry(idx_plot, 7:9), '--')
grid on
hold on
axis equal
%plot(estim(:, 7:9))
plot3(estim(:, 7), estim(:, 8), -estim(:, 9))
title('Measured and Estimated Trajectory')
legend('Measured trajectory','Estimated trajectory')

%attitude
ftime = 0:0.02:(points * 0.02);
figure
plot(exp_telemetry(idx_plot, 10:13),'--')
hold on
grid on
plot(estim(:, 10:13))
title('... Attitude ...')

%linear velocities
figure
plot(ftime', exp_telemetry(idx_plot,1:3),'--')
grid on
hold on
plot(ftime', estim(:,1:3))
title('... Linear Velocities ...')

%angular velocities
size(ftime')
size(exp_telemetry(idx_plot,4:6))
size(estim)
figure
plot(ftime', exp_telemetry(idx_plot,4:6),'--')
grid on
hold on
plot(ftime', estim(:,4:6))
title('... Angular Velocities ...')

%distance from the pivot
figure
dist = sqrt(sum(abs(estim(:,7:9)).^2,2));
dist2 = sqrt(sum(abs(exp_telemetry(idx_plot,7:9)).^2,2));
plot(dist)
grid on
hold on
plot(dist2)
title('... Tether Elongation ...')

%figure
%plot(res_log');
%grid on

%figure
%plot(delays)
%grid on

%%
% Check estimation by forward integration

FORWARD_INT = 1;
if(FORWARD_INT)
    x0_FI = estim(1,7:13);
    %x0_FI = exp_telemetry(start_sample,7:13);
    log_FI = [x0_FI];
    %dt = mean(dt_log);
    RB_INT = rb_model.INT_RED;
    for i = 1:length(estim)-1
    %for i = start_sample:start_sample + points
        vb = estim(i,1:3);
        wb = estim(i,4:6);
        %vb = exp_telemetry(i,1:3);
        %wb = exp_telemetry(i,4:6);
        out = RB_INT(struct('x0', x0_FI, 'p', [vb';wb']));
        x0_FI = full(out.xf);
        log_FI = [log_FI; x0_FI'];
    end
    
    size(flog)
    size(log_FI)
end

%%
% Plot reconstruction results
figure
plot3(exp_telemetry(idx_plot, 7), exp_telemetry(idx_plot, 8), -exp_telemetry(idx_plot, 9))
%plot(exp_telemetry(idx_plot, 10:13))
grid on
hold on
axis equal
%plot(log_FI(:, 4:7),'--')
plot3(log_FI(:, 1), log_FI(:, 2), -log_FI(:, 3))

figure
plot(exp_telemetry(idx_plot, 10:13),'--')
hold on
grid on
plot(log_FI(:, 4:7))