function[X,Y] = ss_solver(sys, init_cond, time, options)
%   [X,Y] = ss_solver(sys, init_cond, time, options) 
%   Solves system of linear differential equations specified in state space:
%   dx/dt = A * x(t) + B * u(t) 
%   y(t) = C * x(t)
%
%   x(t) - state vector [nx1];
%   u(t) - control vector [scalar];
%   y(t) - vector of observations [mx1];
%
%INPUT:
%-sys - describing structure;
% sys.a - dynamics matrix [nxn];
% sys.b - control matrix [nxm];
% sys.c - observation matrix [mxn];
%-init_cond - vector of initial conditions [nx1];
%-time - could a vector of positive real numbers [t0,t1, ... , tn] or range
% [min,max,step] - [seconds]
%-options - is a structure with following fields:
% -- 'time' = 'vector/range' choose between vector format or range.
% -- 'input_type' = 'step/impulse' choose between step responce and impulse
% response
% -- 'control_type' = 'handler/vector/const' could be handler to the input signal function 
% -- 'fundamental_matrix_form' = 'exponential/modal/linear'
% -- 'state_observer' = 'none/luenberger/kalman'
% -- 'process_noise_covariance' [nxn] cov matrix
% -- 'measuring_noise_covariance' [mxm] cov matrix
% -- 'initial_estimation_covariance' [nxn] cov matrix
%
%OUTPUT
%X - vector of solutions
%Y - vector of observations

narginchk(3,5);

%define default options for solver
def_options.time = 'vector';
def_options.input_type = 'step';
def_options.control_type = 'const';
def_options.fundamental_matrix_form = 'modal';
def_options.state_observer = 'none';
def_options.process_noise_covariance = 0;
def_options.measuring_noise_covariance = 0;
def_options.initial_estimation_covariance = 0;
%for more advanced versions of this function
def_options.control = 1;

%check dimensions
init_cond = reshape(init_cond,length(init_cond),1);
if( (size(sys.a,1) ~= size(init_cond,1))  )
    fprintf(2,'sys.a: %d x %d \n', size(sys.a,1), size(sys.a,2));
    fprintf(2, 'initial conditions vector: %d x %d \n', size(init_cond,1), size(init_cond,2));
    error('Initial conditions vector and Dynamics matrix dimesions must agree.');
end

if(size(sys.a,1) ~= size(sys.a,2))
    error('Dynamics matrix should be squared');
end

if(size(sys.b,1) ~= size(sys.a,1))
    fprintf(2, 'Control matrix (b) must have %d rows \n', size(sys.a,1));
    error('Invalid (b) matrix');
end

if(size(sys.c,2) ~= size(sys.a,1))
    fprintf(2, '(c)matrix size [%d x %d] (a)matrix size [%d x %d] \n', size(sys.c,1), size(sys.c,2), size(sys.a,1), size(sys.a,2));
    error('(c) and (a) matrix dimesions must agree');
end

%set solver options
opt = struct;
if(nargin < 4)
   print_options(def_options);
   opt = def_options;
else
    opt = def_options;
    opt_names = fieldnames(options);
    for i=1:length(opt_names)
        prop_name = char(opt_names(i));
        if(isfield(opt, prop_name))
            opt = setfield(opt, prop_name, getfield(options, prop_name));
        end
    end
    print_options(opt);
end

if(isequal(opt.control_type, 'vector') && isequal(opt.time, 'vector'))
    if(~isequal(length(opt.control), length(time)))
        warning('on','Timeline and control vector dimesions do not agree. Output length is set to the shortest of two.');
    end
end

%find a solution
t = []; %time
if (isequal(opt.time, 'vector'))
    t = time;
else if (isequal(opt.time, 'range'))
        if(length(time) > 3)
            error('You chose "range" as a time option, format [min,max,range]');
        end
        t = (time(1):time(3):time(2))';
    else
        error('Unknow option "%s" for the field "time" [vector/range] \n', options.time);
    end
end

%preanalysis
[M,S,stable] = isstable(sys);
I = eye(size(sys.a));
s = diag(S);

%find solution matrix
X = [];
Y = [];
for i=1:length(t)
    %fundamental matrix
    F = eye(size(sys.a));
    if(isequal(opt.fundamental_matrix_form,'exponential'))
        disp('Exponential temporary does not work');
        %stupid crap
        %  |
        % \|/
        F = exp(sys.a * t(i));      
    else if(isequal(opt.fundamental_matrix_form,'modal'))
            if(s == zeros(size(s)))
                F = I;
            else
                F = M * diag( exp(s*t(i)) ) * inv(M);
            end
        else if(isequal(opt.fundamental_matrix_form,'linear'))
                F = fund(sys.a, t(i), 1);
            end
        end
    end
    
    %control
    K = opt.control;
    u = zeros(size(sys.b));
    if(isequal(opt.input_type,'step'))
        if(~isequal(sys.b, zeros(size(sys.b))) || ...
                ~isequal(opt.control,zeros(size(opt.control))) )
            if(isinf(cond(sys.a)))
                if(s == zeros(size(s)))
                    fun = @(x)((I) * sys.b * K);
                else
                    fun = @(x)( (M * diag( exp(s*(t(i) - x)) ) * inv(M))  * sys.b * K);
                end
                int = integral(fun,t(1),t(i),'ArrayValued', true);
                u = F * int;
            else
                u = inv(sys.a) * (F - I) * sys.b * K;
            end
        end
    else
        error('Currently only "step" responce implemented. sorry.');
    end
    
    % solution + observation
    x = F * init_cond + u;
    y = sys.c * x;
    
    X = [X; x'];
    Y = [Y; y'];
end

fprintf('\n                            Simulation done.   \n');


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [] = print_options(options)
   disp('-------------------------');
   disp('ss_solver: Using options:');
   disp('-------------------------');
   disp(options);
end

function [M,S,stable] = isstable(sys)
%returns 
%S,M - eigenvalues and corresponding eigenvectors
%is_stable - stability sign
 
[M,S] = eig(sys.a);
s = diag(S);

eps = 1e-10; %small positive number

positive_poles = 0;
zero_poles = 0;
multiple_poles = 0;
stable = 1;

for i = 1:length(s)
    if(real(s(i)) > 0)
        positive_poles = positive_poles + 1;
    else if (abs(real(s(i))) < eps)
            zero_poles = zero_poles + 1;
        end
    end
    
    if(sum(ismember(s,s(i))) > 1)
        multiple_poles = multiple_poles + 1;
    end
end

if(positive_poles)
    fprintf(2,'\n Warning: \n The system has %d poles with positive real part. \n ', positive_poles);
    stable = 0;
else if(zero_poles)
        fprintf(2,'\n Warning: \n The system has %d poles with zero real parts. \n ', zero_poles);
        stable = 1;
    end
end

if(multiple_poles)
    fprintf(2,'\n Warning: \n The system has %d multiple poles. \n ', multiple_poles);
    stable = 1;
end

%determine system matrix conditioning
large_number = 10000;
if(isinf(cond(sys.a)))
    fprintf(2,'Warning: \n The system matrix is singular to working precision. \n ');
else if(cond(sys.a) > large_number)
         fprintf(2,'Warning: \n The system matrix is ill conditioned \n ');
    end
end

end





