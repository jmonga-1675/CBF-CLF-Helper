%% Example for multiple CBF constraints.
clear all;
close all;
dt = 0.02;
sim_t = 20;
x0 = [0;5;0];

params.v = 1; % velocity
params.u_max = 4; % max yaw rate (left)
params.u_min = -4; % min yaw rate (right)

% Obstacle position
params.xo = [5, 11];
params.yo = [4, 4];
% Obstacle radius
params.d = [2, 2];

% params.xo = 5;
% params.yo = 4;
% params.d = 2;

params.cbf_gamma0 = 1;
% Desired target point
params.xd = 12;
params.yd = 0;

params.clf.rate = 0.5;
% params.weight.slack = 1;
% params.weight.slack = [1, 100, 100];
params.weight.slack = [100, 100];

params.cbf.rate = 50;

dubins = DubinsCar(params);

odeFun = @dubins.dynamics;
% controller = @dubins.ctrlCbfClfQp;
controller = @dubins.ctrlCbfQp;
odeSolver = @ode45;

total_k = ceil(sim_t / dt);
x = x0;
t = 0;   
% initialize traces.
xs = zeros(total_k, dubins.xdim);
ts = zeros(total_k, 1);
us = zeros(total_k-1, 1);
Vs = zeros(total_k-1, 1);
hs = zeros(total_k-1, length(params.xo));
xs(1, :) = x0';
ts(1) = t;
for k = 1:total_k-1
%     t
    % Determine control input.
    % dV_hat: analytic Vdot based on model.
    [u, slack, hs_k, V] = controller(x);        
    us(k, :) = u';
    hs(k, :) = hs_k;
    Vs(k) = V;

    % Run one time step propagation.
    [ts_temp, xs_temp] = odeSolver(@(t, s) odeFun(t, s, u), [t t+dt], x);
    x = xs_temp(end, :)';

    xs(k+1, :) = x';
    ts(k+1) = ts_temp(end);
    t = t + dt;
end

plot_results(ts, xs, us, hs, [params.xo;params.yo], params.d)

function plot_results(t, xs, us, hs, p_o, r_o)

figure
subplot(3,1,1)
plot(t, xs(:,1))
xlabel('t')
ylabel('x [m]')

subplot(3,1,2)
plot(t, xs(:,2))
xlabel('t')
ylabel('y [m]')

subplot(3,1,3)
plot(t, xs(:,3))
xlabel('t')
ylabel('theta [rad]')


figure
plot(t(1:end-1), us)
xlabel('t')
ylabel('u [rad/s]')


lim_min = min(min(xs(:, 1)), min(xs(:, 2)));
lim_max = max(max(xs(:, 1)), max(xs(:, 2)));
lim_min = min([lim_min, p_o(1)-r_o, p_o(2)-r_o]);
lim_max = max([lim_max, p_o(1)+r_o, p_o(2)+r_o]);

figure
plot(xs(:, 1), xs(:, 2)); hold on;
for i = 1:size(p_o, 2)
    draw_circle(p_o(:, i), r_o(i));
end
xlim([lim_min, lim_max]);
ylim([lim_min, lim_max]);
xlabel('x [m]')
ylabel('y [m]')
axis equal

figure
for i = 1:size(hs, 2)
    subplot(size(hs, 2), 1, i)
    plot(t(1:end-1), hs(:, i))
    xlabel('t')
    ylabel('cbf h(s)');
end
end

function h = draw_circle(center,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + center(1);
yunit = r * sin(th) + center(2);
h = plot(xunit, yunit);
hold off

end