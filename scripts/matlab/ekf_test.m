%Extended Kalman Filter test
V = diag([0.01;0.01;0.01; 0.0001;0.005;0.005;0.005].^2);
e_v = [0.5; 0.5; 0.5];             %[0.05; 0.5; 0.1];
e_w = [0.5; 0.5; 0.5];             %[0.1 0.1 0.5]
e_r = [0.5; 0.1; 0.1];             %[0.25; 0.25; 0.05];
e_q = [0.1; 0.05; 0.05; 0.05];    %[0.3; 0.05; 0.05; 0.06];
P_est = 10 * diag([e_v; e_w; e_r; e_q].^2);                 % estimation covariance
W = diag([e_v; e_w; e_r; e_q].^2);

%init condition
x_est = [ -0.3187 0.1826 0.2548 1.8543 -0.1429 -0.1684 -0.2294 -0.0500 -0.7468 0.1894 -0.8363 -0.4818 0.1804];
measurement = [-0.235152, -0.0701073, -0.763765, 0.19925, -0.822125, -0.49634, 0.195082];

rb_model = rigid_body(dt);
rb_model.num.DRE = dre(13, dt);

for k = 1:1
    dt = 0.02;
    %get corresponding control
    u = zeros(3,1);
    %get observation
    z_m = measurement;

    %estimate
    [x_est, P_est, W_est, res] = kiteEKF(x_est, P_est, u, rb_model, dt, W, V, z_m', 'cvodes');
end