% function [xs, us, mus, ts, Vs, dVs_error, etas_zs, feas, slacks, extras] = rollout_controller_eval_clf_FL( ...
%     x0, t0, plant_sys, control_sys, controller, dt, sim_t, with_slack, mu0, verbose, event_options)
function [xs, us, ts, extras] = rollout_controller_eval_clf_FL( ...
    x0, t0, plant_sys, control_sys, controller, dt, varargin)
%% Parse varargin into settings
% Support following symbols for varargin
    % with_slack
    % mu0 = reference for mu input
    % verbose = 0(run quitely), 1(show QP info)
    % event_options = run simulation until event happens
    % sim_t = simulation time
    
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

if ~isfield(settings, 'event_options')
    event_options = 0;
else
    event_options = settings.event_options;
end

if ~isfield(settings, 'sim_t')
    % the dynamics that engage reset map should not be judged by sim_t
    % set sim_t to enough time
    sim_t = 100;
else
    sim_t = settings.sim_t;
end

%% Prepare record paper
% Total time step
total_k = ceil(sim_t / dt);

% Initialize traces.
xs = zeros(1, plant_sys.xdim); % trace of the simulation. 
ys = zeros(1, control_sys.ydim); % Wonsuhk Revised
ts = zeros(2, 1);
Vs = zeros(2, 1);
dys = zeros(1, control_sys.ydim); % Wonsuhk Revised
us = zeros(1, control_sys.udim);
mus = zeros(1, control_sys.udim);
dVs_error = zeros(1, 1);
feas = false(1, 1); % added
slacks = zeros(1, 1); % added
F_sts = zeros(1, 2); % added

x = x0; xs(1, :) = x0;
t = t0; ts(1) = t;

[y, dy, ~, ~, ~] = control_sys.eval_y(x);

z = control_sys.z(x);
ys(1, :) = y; 
dys(1, :) = dy;
Vs(1) = control_sys.clf_FL(y, dy);
zs(1, :) = z;
%% Run Simulation

for k = 1:total_k-1
    % Evaluate CLf related values
    z = control_sys.z(x);
    V = control_sys.clf_FL(y, dy);
    LfV = control_sys.lF_clf_FL(y, dy);
    LgV = control_sys.lG_clf_FL(y, dy);
    
    %% Determine control input
    if k ==1 && ~isempty(mu0)
        % Apply predetermined control input for the initial step.
        mu = mu0; %mu0
        u = control_sys.ctrlFeedbackLinearize(x, mu);
        dV_hat = LfV + LgV * mu; % model based estimate of dV
        feas_k = true;
        slack = 0;
    else
        % Determin control input from the controller.
        [mu, slack, V, feas_k] = controller(x, zeros(control_sys.udim, 1), with_slack, verbose); % mu
        u = control_sys.ctrlFeedbackLinearize(x, mu); % u=u_star + LgLf\mu
        dV_hat = LfV + LgV * mu; % model based estimate of dV
    end
    feas(k, :) = feas_k;
    us(k, :) = u;
    mus(k, :) = mu;
    slacks(k, :) = slack;
    
    % Run simulation for one time step with determined control input.
    if isempty(event_options)
        [ts_k, xs_k] = ode45(@(t, s) plant_sys.dynamics(t, s, u), [t t+dt], x); %TODO => eventfunc!!
    else
        [ts_k, xs_k, te, xe, ~] = ...
            ode45(@(t, s) plant_sys.dynamics(t, s, u), [t t+dt], x, event_options); %TODO => eventfunc!!
    end
    
    t = ts_k(end);
    ts(k+1) = t;
    x_next = xs_k(end, :)';
    % x_next = clip_theta(xs_k(end, :))';
    xs(k+1, :) = x_next;
    
    %% Recording Section: Evaluate CLF related values after one time step.
    % Record dVs_error(k) = dV - dV_hat
    %        Vs = V_next
    [y_next, dy_next, ~, ~, phase_next] = control_sys.eval_y(x_next);
    z_next = control_sys.z(x_next);
    ys(k+1, :) = y_next; % Wonsuhk_revised
    dys(k+1, :) = dy_next;
    V_next = control_sys.clf_FL(y_next, dy_next);
    LfV_next = control_sys.lF_clf_FL(y_next, dy_next);
    LgV_next = control_sys.lG_clf_FL(y_next, dy_next);
    % Estimating Vdot based on second order numerical differentiation.
    dV_hat_next = LfV_next + LgV_next * mu;
    dV_hat = (dV_hat + dV_hat_next)/2;
    % Estimation of actual Vdot based on numerical differentiation
    dV = (V_next - V) / dt;
    dVs_error(k) = dV - dV_hat;
    Vs(k+1) = V_next;
    zs(k+1, :) = z_next;
    
    [F_st_u_next, F_st_nu_next] = plant_sys.get_force(x_next);
    F_st_next = F_st_u_next*u + F_st_nu_next;
    F_sts(k, :) = F_st_next';
        
    if verbose
        fprintf("dV: %.4f \t dV_hat: %.4f \t dV-dV_hat: %.4f \t LfV: %.4f \t theta: %.4f \t dtheta: %.4f\n", ...
            [dV, dV_hat, dV-dV_hat, LfV, x(3), x(4)]);
    end   
    
    if ~isempty(te), break; end
    
    x = x_next;
    y = y_next;
    dy = dy_next;
    % phase = phase_next;
    
end % end of main for loop

% Cache the training data in aggregated form
etas = [ys,dys];
etas_zs = [etas zs];

%% Archive extra datas
% training data required to train the model (refer to GP-CLF-SOCP ver)
extras.training_datas.Vs = Vs;
extras.training_datas.dVs_error = dVs_error;
extras.training_datas.etas_zs = etas_zs;

% Tune this in order to retrieve the specific params you want to save
extras.feas = feas;
extras.slacks = slacks;
extras.mus = mus;
extras.ys = ys;
extras.dys = dys;
extras.Fsts = F_sts;