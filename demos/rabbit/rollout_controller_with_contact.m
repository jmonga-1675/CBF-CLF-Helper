function [xs, us, ts, extras] = rollout_controller_with_contact(...
    x0, plant_sys, control_sys, controller, dt, ...
    event_options, reset_map_function, nstep, ...
    varargin)
%   TODO: to generalize this concept, maybe hybrid simulation might
%   encompass well
%   Input
%       Necessary
%           event_options: function handle that detect the contact
%           reset_map_function: reset the state after detection of the contact
%           nstep: number of contacts
%       Default
%           varargin: with_slack, mu0, verbose, event_options
%
%   Output:
%       xs: trajectory
%       us: control inputs
%       mus: virtual control inputs
%       ts: time stamps
%       etas_zs: [y dy phase dphase]
%       feas: socp feasibility
%       slacks: slack

    settings = parse_function_args(varargin{:});
    
    if ~isfield(settings, 'with_slack')
        with_slack = 0;
    else
        with_slack = settings.with_slack;
    end

    if ~isfield(settings, 'mu0')
        mu0 = [];
    else
        mu0 = settings.mu0;
    end

    if ~isfield(settings, 'verbose')
        verbose = 0;
    else
        verbose = settings.verbose;
    end
    
    %% Record Paper
    x = x0;
    t0 = 0; % TODO: t0 in varargin(?)
    % initialize traces.
    xs = [];
    etas_zs = [];
    ts = [];
    Vs = [];
    us = [];
    mus = [];
    Fsts = [];
    ys = [];
    dys = [];
    dVs_error = [];
    feas = [];
    slacks = [];
    step_logs = [];
    
    %% Run Step Simulation
    for k = 1:nstep
        step_summary = strcat("[Step]", num2str(k)); disp(step_summary);
        
        % Initialize all robot status each step
        is_walking = false;
        
        % change this to lyapunov controller
        [xs_new, us_new, ts_new, extras_new] = rollout_controller_eval_clf_FL( ...
        x0, t0, plant_sys, control_sys, controller, dt, ...
        'with_slack', with_slack, 'verbose', verbose, 'event_options', event_options, 'mu0', mu0);
    
        Vs_new = extras_new.training_datas.Vs;
        dVs_error_new = extras_new.training_datas.dVs_error;
        etas_zs_new = extras_new.training_datas.etas_zs;
    
        ys_new = extras_new.ys;
        dys_new = extras_new.dys;
        Fsts_new = extras_new.Fsts;
        mus_new = extras_new.mus;
        feas_k = extras_new.feas;
        slacks_k = extras_new.slacks;
    
        % TODO: Reset Map when event is detected, no uncertainty term related
        x0 = reset_map_function(xs_new(end, :));
        t0 = ts_new(end);
        
        % If robot has fallen, display the robot's status and stop
        % collecting datas
        height_threshold = control_sys.params.fall_threshold; % threshold to judge if robot has fallen.
        height_trajectory = xs_new(:, 2);
        last_valid = find(height_trajectory<height_threshold, 1, 'first')-1;
        
        if isempty(last_valid)
            last_valid = length(height_trajectory);
            is_walking = true;
        end
        
        % Case I. Robots fell down
        if ~is_walking
            xs_new = xs_new(1:last_valid, :);
            etas_zs_new = etas_zs_new(1:last_valid, :);
            ts_new = ts_new(1:last_valid, :);
            Vs_new = Vs_new(1:last_valid);
            ys_new = ys_new(1:last_valid);
            dys_new = dys_new(1:last_valid);
            
            mus_new = mus_new(1:last_valid-1, :);
            Fsts_new = Fsts_new(1:last_valid-1, :);
            us_new = us_new(1:last_valid-1, :);
            dVs_error_new = dVs_error_new(1:last_valid-1);
            feas_k = feas_k(1:last_valid-1);
            slacks_k = slacks_k(1:last_valid-1);
        end
        
        % Record Simulation data
        xs_new = 0.5*(xs_new(1:end-1, :) + xs_new(2:end, :));
        etas_zs_new = 0.5*(etas_zs_new(1:end-1, :) + etas_zs_new(2:end, :));
        
        ts = [ts; ts_new(1:end-1, :)];
        Vs = [Vs; Vs_new(1:end-1)];
        ys = [ys; ys_new(1:end-1, :)];
        dys = [dys; dys_new(1:end-1, :)];
        
        % training smoothing
        xs = [xs; xs_new];
        etas_zs = [etas_zs; etas_zs_new];
        
        mus = [mus; mus_new];
        us = [us; us_new];
        dVs_error = [dVs_error; dVs_error_new'];
        feas = [feas; feas_k];
        slacks = [slacks; slacks_k];
        Fsts = [Fsts; Fsts_new];
        step_logs = [step_logs; size(mus, 1)];

        if ~is_walking
            disp("Robot fell down!");
            break;
        end
    end
    
    extras.training_datas.Vs = Vs;
    extras.training_datas.dVs_error = dVs_error;
    extras.training_datas.etas_zs = etas_zs;
    
    extras.feas = feas;
    extras.slacks = slacks;
    extras.mus = mus;
    extras.ys = ys;
    extras.dys = dys;
    extras.Fsts = Fsts;
    extras.Fsts = Fsts;
    extras.step_logs = step_logs;
end

