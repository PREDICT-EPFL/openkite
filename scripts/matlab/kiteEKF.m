function[x_est, P_est, W_est, res] = kiteEKF(x, P, u, model, dt, W, V, z_m, method)
%EXTENDED Kalman filter for Kite state estimation

%Measurement matrix
H = [zeros(7,6) eye(7,7)];

%propagation step
%state estimate
if(isequal(method, 'rk4'))
    x_est = model.num.INT({x, u, dt});
    x_est = full(x_est{1});
else if (isequal(method, 'cvodes'))
        %out = model.num.INT(struct('x0', x, 'p',u));
        out = model.INT(struct('x0', x));
        x_est = full(out.xf);       
    end
end

%covariance estimate
I = eye(length(x_est));

%A = model.num.DYN_JACOBIAN({x, u});
%simple rigid body model
A = model.JACOBIAN({x, u});
%compute fundamental matrix
%[EM, ES] = eig( full(A{1}) );
%es = diag(ES);

%F = EM * diag(exp(es * dt)) / EM;
A = full(A{1});
F = eye(13) + A * dt;
%A = model.num.RK4_JACOBIAN({x, u, dt});
%F = full(A{1});
%F = x_est * (1/x);
P_tmp = (F * P) * F';
P_est = (F * P) * F' + W;

%use RK4 integration routine
P_est_vec = model.num.DRE(struct('x0', vec(P), 'p', [vec(A); vec(W)]));
P_est = reshape(full(P_est_vec.xf), size(A));
norm(P_est,2);

%adaptation
gamma = 0.01;
x_pred = x_est;
P_pred = P_est;

%observability test
%modes = eig(F);
%O = [H', F' - diag(modes)];
%if (rank(O) < 13)
%    warning('Low observability')
%end

%update step
y = z_m - H * x_est;
res = y;
if(norm(y,2) >= 0.25)
    W_est = W;
    return
end
%innovation covariance
S = (H * P_est) * H' + V;
%Kalman gain
K = (P_est * H') / S;

%updated state estimate
x_est = x_est + K * y;
x_est = real(x_est);
%estimate covariance update
P_est = (I - K * H) * P_est;

%noise adaptation step
r = x_est - x_pred;
Q_star = r * r' + P_pred - P_est - W;
W_est = gamma * Q_star + (1 - gamma) * W;
%[E,S] = eig(W_est);
%W_est = real(S);
end




