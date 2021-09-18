%% RABBIT SIMULATION
% Wonsuhk Jung smw04015@snu.ac.kr
% This is CBF-CLF Helper Version of RABBIT Bipedal Walker.
% You can use a controller designed based on RES-CLF
close all; clear all; clc;

%% Animation setting
% This is only designed for enabling animation of bipedal walker
global animation_scale 

%% Set Parameters
% Simulation Setting
    % sim_t: overall simulation time: do not need it in this mode actually
    % dt: time interval to collect data used for plot/animation
    % nstep: number of step that RABBIT would take
% params: all the system-related parameters
% verbose: speak about the CBF-CLF(1), silent(0)
% with slack: slack variable in QP(1), no slack variable in QP(0)

params = init_clf_simulation_rabbit; % Load all the parameters
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

%% Prepare System (Uncertainty control)
control_sys = RabbitBuiltIn(params);
plant_sys = RabbitBuiltIn(params);

% Reflect model uncertainty here
plant_sys.params.scale = 1.5;
%plant_sys.params.torso_add = 10;

controller = @control_sys.ctrlClfQpFL;
event_options = odeset('Events', @plant_sys.rabbit_event);

%% Prepare record paper
xs = [];
ts = [];
us = [];
extras = [];

% Initialize the state and time
x = x0;
t0 = 0;
%% Main Simulation
for k = 1:nstep
    step_summary = strcat("[Step]", num2str(k)); disp(step_summary);
    
    [xs_new, us_new, ts_new, extras_new] = rollout_controller_eval_clf_FL( ...
    x0, t0, plant_sys, control_sys, controller, dt, ...,
    'event_options', event_options, 'sim_t', sim_t, 'with_slack', with_slack, ...,
    'verbose', verbose);
    
    % Reset Map when event is detected
    n = plant_sys.params.xdim/2;
    q = R * xs_new(end, 1:n)';
    dq = R*plant_sys.dq_pos_gen_v2([xs_new(end, 1:n) xs_new(end, n+1:2*n)]);
    
    % Renew the new state
    x0 = [q; dq];
    t0 = ts_new(end);
    
    ts = [ts; ts_new(1:end-1)];
    xs = [xs; xs_new(1:end-1, :)];
    us = [us; us_new];
    extras = [extras; extras_new];
end
%% Plot the result
% TODO: tune the plot setting into new version 
plot_rabbit_result(xs, ts, us, extras);

%% Five Link Animation
animation_scale = plant_sys.params.scale;
animation_dt = 0.05;
x_quan_vec = coordinateTransformation(xs);
anim_flat_ground(ts, x_quan_vec, animation_dt)