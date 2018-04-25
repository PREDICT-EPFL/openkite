function [center, radius, rmse] = sphere_fit_3d(points, icenter, iradius)
max_iter = 500;
gamma = 2 * 1e-6;
x = [icenter;iradius];
for i = 1:max_iter
    cost = cost(x, points);
    gradient = gradient(x, points);
    x = x - gamma * gradient;
    fprintf('[%f] [%f %f %f %f] [%f %f %f %f] \n', cost, gradient(1), ...
        gradient(2), gradient(3), gradient(4), x(1), x(2), x(3), x(4));
end
center = x(1:3);
radius = x(4);
rmse = sqrt(cost);
end

function [cost] = cost(param, points)
cost = 0;
center = param(1:3);
radius = param(4);
for i = 1:size(points, 1)
    cost = cost + ( norm(points(i, :)' - center, 2)^2 - radius^2 )^2;
end
cost = 0.5 * cost;
end

function [gradient] = gradient(param, points)
gradient = [0;0;0;0];
center = param(1:3);
radius = param(4);

for i = 1:size(points, 1)
    Li = norm(center - points(i,:)', 2)^2 - radius^2;
    gradient = gradient - 2 * Li * [points(i ,:)' - center; radius];
end
end