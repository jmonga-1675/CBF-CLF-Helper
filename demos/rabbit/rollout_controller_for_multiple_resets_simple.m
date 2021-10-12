function [xs, us, ts, extras] = rollout_controller_for_multiple_resets_simple(...
    x0, plant_sys, control_sys, controller, dt, ...
    event_func, reset_map_function, nstep, ...
    varargin)
%   Rollout control for the hybrid system that uses reset_map functions
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
    
    if ~isfield(settings, 'mu_ref')
        mu_ref = [];
    else
        mu_ref = settings.mu_ref;
    end

    if ~isfield(settings, 'verbose')
        verbose = 0;
    else
        verbose = settings.verbose;
    end
    
    %% Record Paper
    sim_t = inf; % set enough time if sim_t is not the standard of the end
    x = x0;
    t0 = 0; % TODO: t0 in varargin(?)
    % initialize traces.
    xs = [];
    ts = [];
    Vs = [];
    us = [];
    mus = [];
    feas = [];
    slacks = [];
    
%     FL_linearizer = @(obj, x, mu) (obj.ctrlFeedbackLinearize(x, mu));
%     FL_logger = @(obj, x, mu) (obj.logRecords(x, mu));
    %% Run Step Simulation
    for k = 1:nstep
        step_summary = strcat("[Step]", num2str(k)); disp(step_summary);
        
        % Initialize all robot status each step
        is_walking = false;
        
        [xs_new, us_new, ts_new, extras_new] = rollout_controller( ...
        x0, plant_sys, control_sys, controller, sim_t,...
        't0', t0, 'dt', dt, 'end_event', event_func, ...,
        'verbose', verbose, 'with_slack', with_slack);
    
        Vs_new = extras_new.Vs;
        mus_new = extras_new.mus;
        feas_k = extras_new.feas;
        slacks_k = extras_new.slacks;
    
        % TODO: Reset Map when event is detected, no uncertainty term related
        x0 = reset_map_function(xs_new(:, end)');
        t0 = ts_new(end);
        
        % If robot has fallen, display the robot's status and stop
        % collecting datas
        height_threshold = control_sys.params.fall_threshold; % threshold to judge if robot has fallen.
        height_trajectory = xs_new(2, :);
        last_valid = find(height_trajectory<height_threshold, 1, 'first')-1;
        
        if isempty(last_valid)
            last_valid = length(height_trajectory);
            is_walking = true;
        end
        
        % Case I. Robots fell down
        if ~is_walking
            xs_new = xs_new(:, 1:last_valid); 
            ts_new = ts_new(:, 1:last_valid);
            Vs_new = Vs_new(:, 1:last_valid);
            
            mus_new = mus_new(:, 1:last_valid-1);
            us_new = us_new(:, 1:last_valid-1);
            feas_k = feas_k(1:last_valid-1);
            slacks_k = slacks_k(1:last_valid-1);
        end
        
        % Record Simulation data
        xs_new = 0.5*(xs_new(:, 1:end-1) + xs_new(:, 2:end));
        
        ts = [ts, ts_new(:, 1:end-1)];
        Vs = [Vs, Vs_new(1:end-1)];
        
        % training smoothing
        xs = [xs, xs_new];
        mus = [mus, mus_new];
        us = [us, us_new];
        feas = [feas, feas_k];
        slacks = [slacks, slacks_k];

        if ~is_walking
            disp("Robot fell down!");
            break;
        end
    end
    
    extras.logs.Vs = Vs;
    extras.feas = feas;
    extras.slacks = slacks;
    extras.mus = mus;
end

