%% RABBIT SIMULATION
% Wonsuhk Jung smw04015@snu.ac.kr
% This is CBF-CLF Helper Version of RABBIT Bipedal Walker.
% You can use a controller designed based on RES-CLF
close all; clear all; clc;

%% Animation setting
% This is only designed for enabling animation of bipedal walker
global animation_scale 

%% Step 1. Set Simulation Setting and System-parameters
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

%% Step 2. Prepare System (Uncertainty control)
control_sys = RabbitBuiltIn(params);
plant_sys = RabbitBuiltIn(params);

%% Reflect model uncertainty here
plant_sys.params.scale = 1.5;
%plant_sys.params.torso_add = 10;

% controller = @control_sys.ctrlClfQpFL; % controller
feedback_clf_qp_controller = @(x, mu_ref, with_slack, verbose) control_sys...
        .wrap_controller(x, mu_ref, with_slack, verbose, @control_sys.ctrlClfQpFL);
simple_controller = feedback_clf_qp_controller;
complex_controller = @control_sys.ctrlClfQpFL;

event_func = @plant_sys.rabbit_event; % contact detector
reset_map_function = @plant_sys.reset_map; % reset mapper

%% Step 3. Initialize state
% Init state.
q0_ini = params.legacy.q0_ini;
dq0_ini = params.legacy.dq0_ini;
x0 = [q0_ini;dq0_ini];
x = x0;
t0 = 0;

%% Step 4. Main Simulation
% rollout_controller_with_contact_version(for hybrid system)

% run_simple
%   to exploit full information, set this to false
%   if you only want to use mu, y, set this to true

% TODO: rollout_controller_for_multiple_resets.m => transpose different
    % use eval_FL
run_simple = true;
if run_simple
    [xs, us, ts, extras] = rollout_controller_for_multiple_resets_simple(...
        x0, plant_sys, control_sys, simple_controller, dt, ...
        event_func, reset_map_function, nstep,...
        'with_slack', with_slack, 'verbose', verbose);
    xs = xs';
    ts = ts';
else
    [xs, us, ts, extras] = rollout_controller_for_multiple_resets_complex(...
        x0, plant_sys, control_sys, complex_controller, dt, ...
        event_func, reset_map_function, nstep,...
        'with_slack', with_slack, 'verbose', verbose);
end

%% Step 5. Plot the result
if ~run_simple
    plot_rabbit_result(xs, ts, us, extras);
end

%% Step 6. Five Link Animation
animation_scale = plant_sys.params.scale;
animation_dt = 0.05;
x_quan_vec = coordinateTransformation(xs);
anim_flat_ground(ts, x_quan_vec, animation_dt)