function [xs, us, ts, Vs, dVs_error, ys] = rollout_feedback_linearization( ...
    x0, plant_sys, control_sys, controller, dt, sim_t, with_slack, mu0, verbose)
if nargin < 7
    with_slack = 0;
elseif isempty(with_slack)
    with_slack = 0;
end

if nargin < 8
    mu0 = [];
end

if nargin < 9
    verbose = 0; % run queitly 
end

% Total time step
total_k = ceil(sim_t / dt);

% Initialize traces.
xs = zeros(total_k, plant_sys.xdim); % trace of the simulation. 
ys = zeros(total_k, 1);
us = zeros(total_k - 1, control_sys.udim);
ts = zeros(total_k, 1);
Vs = zeros(total_k, 1);
dVs_error = zeros(total_k - 1, 1);

x = x0; xs(1, :) = x0;
t = 0; ts(1) = t;

[y, dy, ~, ~, phase] = control_sys.eval_y(x)
ys(1) = y;
Vs(1) = control_sys.clf_FL(y, dy);
%% Run Simulation
for k = 1:total_k-1
    % Evaluate CLf related values
    V = control_sys.clf_FL(y, dy);
    LfV = control_sys.lF_clf_FL(y, dy);
    LgV = control_sys.lG_clf_FL(y, dy);
    
    %% Determine control input
    if k ==1 && ~isempty(mu0)
        % Apply predetermined control input for the initial step.
        mu = mu;
        u = control_sys.ctrlFeedbackLinearize(x, mu);
        dV_hat = LfV + LgV * mu; % model based estimate of dV
    else
        % Determin control input from the controller.
        mu = controller(x, zeros(control_sys.udim, 1), with_slack, verbose);        
        u = control_sys.ctrlFeedbackLinearize(x, mu);
        %% for debug;
%         u = 0;
%         feedforward = -control_sys.lglf_y(x)\control_sys.l2f_y(x);
%         mu = control_sys.lglf_y(x) * (u - feedforward);        
        dV_hat = LfV + LgV * mu; % model based estimate of dV
    end
    us(k, :) = u;
    
    % Run simulation for one time step with determined control input.
    [ts_k, xs_k] = ode45(@(t, s) plant_sys.dynamics(t, s, u), [t t+dt], x);
    t = ts_k(end);
    ts(k+1) = t;
    x_next = xs_k(end, :)';
    % x_next = clip_theta(xs_k(end, :))';
    xs(k+1, :) = x_next;
    
    % Evaluate CLf related values after one time step.
    [y_next, dy_next, ~, ~, phase_next] = control_sys.eval_y(x_next);
    ys(k+1, 1) = y_next;
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
%     if abs(dV-dV_hat) > 50
%         disp("debug");
%     end
    if verbose
        fprintf("mu: %.4f \t dV: %.4f \t dV_hat: %.4f \t dV-dV_hat: %.4f \t LfV: %.4f \t LgV: %.4f \t theta: %.4f \t dtheta: %.4f\n", ...
            [mu, dV, dV_hat, dV-dV_hat, LfV, LgV, x(3), x(4)]);
    end   
    
    x = x_next;
    y = y_next;
    dy = dy_next;
    phase = phase_next;
end % end of main for loop