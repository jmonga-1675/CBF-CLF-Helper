%% RABBIT SIMULATION
% This is CBF-CLF Helper Version of RABBIT Bipedal Walking.
% This version now supports
close all; clear all; clc;
%% Set Parameters
params = init_clf_simulation_rabbit; % save

% Simulation Setting
    % sim_t: overall simulation time: do not need it in this mode actually
    % dt: time interval to collect data used for plot/animation
    % nstep: number of step that RABBIT would take
dt = 0.025; % simulation interval
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
controller = @control_sys.ctrlClfQpFL;
event_func = @rabbit_event;
%odeFun = @plant_sys.dynamics;
%odeSolver = @ode45;
%options = odeset('Events', @rabbit_event);

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
u_prev = [0;0]; % Todo

%% Main Simulation
for k = 1:nstep
    step_summary = strcat("[Step]", num2str(k));
    disp(step_summary);
    [xs_new, us_new, ts_new, Vs_new, ~, ys_new] = rollout_feedback_linearization( ...
    x0, t0, plant_sys, control_sys, controller, dt, sim_t, with_slack, [], verbose, event_func);
    
    % Reset Map when event is detected
    n = params.xdim/2;
    q = R * xs_new(end, 1:n)';
    dq = R*dq_pos_gen_v2([xs_new(end, 1:n) xs_new(end, n+1:2*n)], plant_sys.params.scale);
    
    % Renew the new state
    x0 = [q; dq];
    t0 = ts_new(end);
    
    ts = [ts; ts_new(2:end)];
    xs = [xs; xs_new(2:end, :)];
end

%% Five Link Animation
animation_dt = 0.05;
x_quan_vec = coordinateTransformation(xs);
anim_flat_ground(ts, x_quan_vec, animation_dt)
