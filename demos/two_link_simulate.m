function [xRec, fRec, uRec, tRec, xRec_tilde] = two_link_simulate(x0, beta, Nsteps, th1d, sim_option, sim_event_option)
%TWO_LINK_SIMULATE simulates the two link robot for Nsteps number of
%walking steps and with I/O parameters beta and initial condition x0.
% Ayush Agrawal, Jason Choi
% ayush.agrawal@berkeley.edu
% jason.choi@berkeley.edu
% Inputs:   x0: initial state of the step in full representation in the original coordinate (8 x 1 vector)
%           beta: bezier coefficients that define the reference gait
%           Nsteps: number of walking steps to simulate
%           th1d: desired value of q1 (phase var) at the end of a step
%           sim_option: simulation option.
%              'vanilla': use two_link_dynamics.m to run simulation.
%               'tilde': use two_link_dynamics_tilde.m to run simulation.
%               'vanilla_class': use CompassWalker.dynamics to run simulation.
%               'tilde_class': use CompassWalkerTilde.dynamics to run simulation.
%           sim_event_option: event used for reset map in simulation.
%               'angle': Reset is defined as q1 reaching q1_desired.
%               'foot': Reset is defined as foot position reaching the
%               ground (y=0).

if nargin<1 || isempty(x0)
    x0 = getInitialCondition();
end
if nargin<3
    Nsteps = 1;
end
if nargin<4
    th1d = -x0(3);
end
if nargin<5
    sim_option = 'vanilla';
end
if nargin<6
    sim_event_option = 'foot';
end

% variables to record data
tRec = 0;
xRec = [];
xRec_tilde = [];
% parameters
params.ctrl_type = 'io';
params.beta = beta;
params.eps = 0.1; % Epsilon used in feedback linearization.
params.th1d = th1d;
params.reset_q1_threshold = -0.03;

% Reflection matrix (swap leg1 and 2)
R = two_link_reflect;

for i = 1:Nsteps
    step_summary.x_init = x0;
    [xs, ts, xs_tilde] = step_simulate_helper(x0, tRec(end), params, sim_option, sim_event_option);    
    
    % Record data
    tRec = [tRec;ts];
    xRec = [xRec;xs];
    xRec_tilde = [xRec_tilde;xs_tilde];
    step_summary.x_terminal = xs(end, :)';
    % Check swing foot position
    pSw_at_event = pSw_gen(xs(end, :)');
    step_summary.pSw_terminal = pSw_at_event;
    
    print_step_summary(i, step_summary);
    
    %% Apply reset map and reinitialize initial state for next step.
    q_minus = xs(end, 1:4)';
    dq_minus = xs(end, 5:end)';
    q_plus = R*q_minus;
    dq_plus = R*dqPlus_gen([q_minus;dq_minus]);
    x0 = [q_plus;dq_plus];    
end

tRec = tRec(2:end);
uRec = [];
fRec = [];
%% Compute control inputs and Contact forces
for i = 1:length(tRec)
    u = two_link_io_control(xRec(i, :)', params);
    uRec = [uRec; u];
    fRec = [fRec; Fst_gen(xRec(i, :)', u)'];
end
end

function [xs, ts, xs_tilde] = step_simulate_helper(x0, t0, params, sim_option, event_option)
%% function that runs walk step simluation
% Inputs:   x0: initial state of the step in full representation in the original coordinate (8 x 1 vector)
%           t0: initial time of the step
%           params: control params
%               .ctrl_type
%               .beta
%               .eps
%               .th1d
%           sim_option: simulation option
%              'vanilla': use two_link_dynamics.m to run simulation.
%              'tilde': use two_link_dynamics_tilde.m to run simulation.
%              'vanilla_class': use CompassWalker.dynamics to run simulation.
%              'tilde_class': use CompassWalkerTilde.dynamics to run simulation.
    tspan = [t0, t0+3];
    gait_params = [params.beta;params.th1d];
    if strcmp(sim_option, 'vanilla')
        dynamics_function = @(t,x)two_link_dynamics(t, x, params);
        event_function = @(t, x)two_link_event(t, x, params, event_option);
        odeOpts = odeset('Events', event_function);
        [ts, xs] = ode45(dynamics_function, tspan, x0, odeOpts);
        xs_tilde = xs_tilde_from_original(xs, gait_params);           
    elseif strcmp(sim_option, 'vanilla_class')
        dynsys = CompassWalker(1:4, gait_params, [], event_option);
        dynamics_function = @(t,x) dynsys.dynamics(t, x, [], [], true, params);
        event_function = @dynsys.reset_event;
        odeOpts = odeset('Events', event_function);
        x0_condensed = x0([3, 4, 7, 8]);
        init_stance_position = pSt_gen(x0);
        [ts, xs_condensed] = ode45(dynamics_function, tspan, x0_condensed, odeOpts);
        xs = convert_condensed_states_to_full_states( …
            xs_condensed, init_stance_position);
        xs_tilde = xs_tilde_from_original(xs, gait_params);
    elseif strcmp(sim_option, 'tilde')
        x0_tilde = [q_tilde_from_q_gen(x0(1:4), gait_params); 
            dq_tilde_from_qdq_gen(x0(1:4), x0(5:8), gait_params)];
        dynamics_function = @(t,x)two_link_dynamics_tilde(t, x, params);
        event_function = @(t, x)two_link_event_tilde(t, x, params, event_option);
        odeOpts = odeset('Events', event_function);
        [ts, xs_tilde] = ode45(dynamics_function, tspan, x0_tilde, odeOpts);
        xs = xs_original_from_tilde(xs_tilde, gait_params);
    elseif strcmp(sim_option, 'tilde_class')
        dynsys = CompassWalkerTilde(1:4, gait_params, [], event_option);
        dynamics_function = @(t,x) dynsys.dynamics(t, x, [], [], true, params);
        event_function = @dynsys.reset_event;
        odeOpts = odeset('Events', event_function);
        init_stance_position = pSt_gen(x0);
        x0_tilde = [q_tilde_from_q_gen(x0(1:4), gait_params); 
            dq_tilde_from_qdq_gen(x0(1:4), x0(5:8), gait_params)];
        x0_condensed_tilde = x0_tilde([3, 4, 7, 8]);
        [ts, xs_condensed_tilde] = ode45( …
            dynamics_function, tspan, x0_condensed_tilde, odeOpts);
        xs_tilde_temp = zeros(length(ts), 8);
        xs_tilde_temp(:, [3,4,7,8]) = xs_condensed_tilde;
        xs_original_temp = xs_original_from_tilde(xs_tilde_temp, gait_params);
        xs_condensed = xs_original_temp(:, [3, 4, 7, 8]);
        xs = convert_condensed_states_to_full_states( …
            xs_condensed, init_stance_position);
        xs_tilde = xs_tilde_from_original(xs, gait_params);
    end
end

function print_step_summary(i, log)
    fprintf("———Summary of %d-th step———\n", i);
    disp("Initial state");
    disp(log.x_init');
    disp("Terminal state");
    disp(log.x_terminal');
    disp("Swing foot position at terminal state (x, y)");
    disp(log.pSw_terminal');    
end