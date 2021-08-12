%% RABBIT SIMULATION
% This is CBF-CLF Helper Version of RABBIT Bipedal Walking.
% This version now supports
close all; clear all; clc;
%% Animation setting
global animation_scale % animation scale declared global

%% Set Parameters
% Simulation Setting
    % sim_t: overall simulation time: do not need it in this mode actually
    % dt: time interval to collect data used for plot/animation
    % nstep: number of step that RABBIT would take
% params: all the system-related parameters
% verbose: speak about the CBF-CLF(1), silent(0)
% with slack: slack variable in QP(1), no slack variable in QP(0)

params = init_clf_simulation_rabbit; % Load all the parameters
dt = 0.01; % simulation interval
sim_t = 100; % simulation time 
nstep = 10; % number of steps
with_slack = true; % use slack variable in QP (change this into false in order not to use it)
verbose = false; % show the QP-solving process

%% Initialize state
% Init state.
q0_ini = params.legacy.q0_ini;
dq0_ini = params.legacy.dq0_ini;
x0 = [q0_ini;dq0_ini];
% Relabeling matrix (reset map)
R = [1, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 0, 1;
     0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0];

%% Prepare System
control_sys = RabbitBuiltIn(params);
plant_sys = RabbitBuiltIn(params);
plant_sys.params.scale = 3; % Reflect model uncertainty here
%plant_sys.params.torso_add = 10;
controller = @control_sys.ctrlClfQpFL;
event_options = odeset('Events', @plant_sys.rabbit_event);
%odeFun = @plant_sys.dynamics;
%odeSolver = @ode45;

%% Prepare record paper
total_k = ceil(sim_t / dt);
x = x0;
t0 = 0;   
% initialize traces.
% xs = zeros(total_k, dynsys.xdim);
xs = [];
ts = [];
us = zeros(total_k-1, control_sys.udim);
hs = zeros(total_k-1, 1);
Vs = zeros(total_k-1, 1);
Fsts = zeros(total_k-1, 2);

xs(1, :) = x0';
ts(1) = t0;

%% Main Simulation
for k = 1:nstep
    step_summary = strcat("[Step]", num2str(k));
    disp(step_summary);
    
    [xs_new, us_new, ts_new, Vs_new, ~, ys_new] = rollout_feedback_linearization( ...
    x0, t0, plant_sys, control_sys, controller, dt, sim_t, with_slack, [], verbose, event_options);
    
    % Reset Map when event is detected
    n = params.xdim/2;
    q = R * xs_new(end, 1:n)';
    dq = R*plant_sys.dq_pos_gen_v2([xs_new(end, 1:n) xs_new(end, n+1:2*n)]);
    
    % Renew the new state
    x0 = [q; dq];
    t0 = ts_new(end);
    
    ts = [ts; ts_new(2:end)];
    xs = [xs; xs_new(2:end, :)];
end
%% New Figures
subplot_gap = [0 0.1];
fig_sz = [8 3]; 
fig_sz_with_subplots = [8 4]; 
plot_pos = [0 0 8 3];

main_color_orange = [217 83 25]/255; % r,g,b
main_color_dark_red = [162 20 47]/255; % r,g,b
main_color_blue = [0 114 189]/255; % r,g,b
main_color_green = [119 172 48]/255; 

figure(1);
plot(ts(2:end), plant_sys.y_out_history(:,1)); grid on;
xlabel('Time(s)'); ylabel('q_1');
title('Stance Leg Hip Angle');

figure(2);
plot(ts(2:end),plant_sys.y_out_history(:,2)); grid on;
xlabel('Time(s)'); ylabel('q_3');
title('Stance Leg Knee Angle');

figure(3);
plot(ts(2:end),plant_sys.y_out_history(:,3)); grid on;
xlabel('Time(s)'); ylabel('q_2');
title('Swing Leg Hip Angle');

figure(4);
plot(ts(2:end),plant_sys.y_out_history(:,4)); grid on;
xlabel('Time(s)'); ylabel('q_4');
title('Swing Leg Knee Angle');

figure(5); plot(ts(2:end), plant_sys.F_st_history(:,1)); grid on; hold on; xlabel('t');
ylabel('FSt'); plot(ts(2:end), plant_sys.F_st_history(:,2)); legend(['FSt(1)';'FSt(2)']);

figure(6);
plot(ts(2:end),abs(plant_sys.F_st_history(:,1)./plant_sys.F_st_history(:,2)), 'color', main_color_blue, 'LineWidth', 1); grid on; hold on;
title('Contact Force Ratio ($|F_{T}$/$F_{N}|$)', 'Interpreter', 'latex');
xlim([0 ts(end)+1]);
xlabel('Time (sec)');
ylabel('$|F_{T}$/$F_{N}|$', 'Interpreter', 'latex');
callPlotSettings(fig_sz, plot_pos);
%print('-f6','plot_F_ratio', '-dpng');

figure(7);
plot(ts(2:end),plant_sys.u_motor_history); grid on; hold on;
title('torques'); legend('u1','u2','u3','u4');

figure(9); 
subplot1(2,1,'Gap', subplot_gap);
subplot1(1); hold on;
title('Euclidian Norm of the Tracking Error');
plot(ts(2:end),sqrt(plant_sys.y_out_history(:,1).^2+plant_sys.y_out_history(:,2).^2+...
    plant_sys.y_out_history(:,3).^2+plant_sys.y_out_history(:,4).^2), 'color', main_color_blue); 
ylabel('||y||_2 (rad)'); grid on;
callPlotSettings(fig_sz_with_subplots, plot_pos);

subplot1(2); hold on;
title('Euclidian Norm of the Derivative of the Tracking Error');
plot(ts(2:end),sqrt(plant_sys.dy_out_history(:,1).^2+plant_sys.dy_out_history(:,2).^2+...
    plant_sys.dy_out_history(:,3).^2+plant_sys.dy_out_history(:,4).^2), 'color', main_color_blue); 
ylabel('$||\dot{y}||_{2}$ (rad/s)','Interpreter','latex'); xlabel('Time (s)'); grid on;
callPlotSettings(fig_sz_with_subplots, plot_pos);
%print('-f9','plot_tracking_error', '-dpng');
%% Five Link Animation
animation_scale = plant_sys.params.scale;
animation_dt = 0.05;
x_quan_vec = coordinateTransformation(xs);
anim_flat_ground(ts, x_quan_vec, animation_dt)
