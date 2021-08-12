function [xs, us, ts, Vs, dVs_error, ys] = rollout_feedback_linearization( ...
    x0, t0, plant_sys, control_sys, controller, dt, sim_t, with_slack, mu0, verbose, event_options)
% controller = CBF-CLF-QP
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

if nargin < 10
    event_options = []; % 0: no event function
end

% Total time step
total_k = ceil(sim_t / dt);

% Initialize traces.
xs = zeros(2, plant_sys.xdim); % trace of the simulation. 
ys = zeros(2, control_sys.ydim); % Wonsuhk Revised
us = zeros(1, control_sys.udim);
ts = zeros(2, 1);
Vs = zeros(2, 1);
dVs_error = zeros(1, 1);

x = x0; xs(1, :) = x0;
t = t0; ts(1) = t; % wonsuhk revised

[y, dy, ~, ~, phase] = control_sys.eval_y(x);
ys(1, :) = y;
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
        mu = mu0; %mu0
        u = control_sys.ctrlFeedbackLinearize(x, mu);
        dV_hat = LfV + LgV * mu; % model based estimate of dV
    else
        % Determin control input from the controller.
        mu = controller(x, zeros(control_sys.udim, 1), with_slack, verbose);        
        u = control_sys.ctrlFeedbackLinearize(x, mu); % u=u_star + LgLf\mu
        %% for debug;
%         u = 0;
%         feedforward = -control_sys.lglf_y(x)\control_sys.l2f_y(x);
%         mu = control_sys.lglf_y(x) * (u - feedforward);        
        dV_hat = LfV + LgV * mu; % model based estimate of dV
    end
    us(k, :) = u;
    
    % Run simulation for one time step with determined control input.
    if isempty(event_options)
        [ts_k, xs_k] = ode45(@(t, s) plant_sys.dynamics(t, s, u), [t t+dt], x); %TODO => eventfunc!!
    else
        [ts_k, xs_k, te, ~, ~] = ode45(@(t, s) plant_sys.dynamics(t, s, u), [t t+dt], x, event_options); %TODO => eventfunc!!
    end
    
    t = ts_k(end);
    ts(k+1) = t;
    x_next = xs_k(end, :)';
    % x_next = clip_theta(xs_k(end, :))';
    xs(k+1, :) = x_next;
    
    % Recording Section: Evaluate CLF related values after one time step.
    % Record dVs_error(k) = dV - dV_hat
    %        Vs = V_next
    [y_next, dy_next, ~, ~, phase_next] = control_sys.eval_y(x_next);
    ys(k+1, :) = y_next; % Wonsuhk_revised
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
        mu'
        fprintf("dV: %.4f \t dV_hat: %.4f \t dV-dV_hat: %.4f \t LfV: %.4f \t theta: %.4f \t dtheta: %.4f\n", ...
            [dV, dV_hat, dV-dV_hat, LfV, x(3), x(4)]);
        LgV
    end   
    
    if ~isempty(te)
        disp("Touched Ground");
        break;
    end
    
    x = x_next;
    y = y_next;
    dy = dy_next;
    phase = phase_next;
end % end of main for loop