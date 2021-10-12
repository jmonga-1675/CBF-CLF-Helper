function [xs, us, ts, extraout] = rollout_controller( ...
    x0, plant_sys, control_sys, controller, T, varargin)
%% [xs, us, ts, extras] = rollout_controller( ...
%%    x0, t0, plant_sys, control_sys, controller, T, varargin)
%% [xs, us, ts, extras] = rollout_controller( ...
%%    x0, t0, plant_sys, control_sys, controller, T, ...
%%    'field1_name', field1_value, 'field2_name', field2_value)
%% example:
%% [xs, us, ts, extraout] = rollout_controller( ...
%%    x0, t0, plant_sys, control_sys, @control_sys.ctrlClfQp, T);
%% Vs = extraout.Vs;
%% simulates plant_sys with ctrlClfQp controller of the control_sys
%% for [t0, t0+T] starting at the initial state x0.
%% INPUTS
%   x0: initial state
%   plant_sys: ControlAffineSys instance for the plant simulation.
%   control_sys: ControlAffineSys instance for the model and the controller.
%   controller: function handle for controller.
%       examples: ctrlClfQp, ctrlCbfQp, ctrlCbfClfQp
%   T: simulation time (time horizon)
%% Supported fields for varargin
%   SIMULATION OPTIONS:
%   t0: initial t (default: 0.0)
%   u0: initial control input (if you want to specify the first input.)
%   dt: sampling time (default: 0.01)
%   end_event: function handle that ends the simulation, for more
%   information, please check out matlab ode Events.
%   ode_func: your own choice of ode_func for simulation. (default: ode45)
%   CONTROLLER OPTIONS:
%   with_slack: if 1: relax constraints, 0: hard constraints. (default 1)
%   ratio_u_diff: ratio(0.0~1.0) in the min-norm objective for minimizing the
%   difference between u and u_prev (default 0.0)
%       if 0.0 : minimizes the norm of plain u.
%       if 1.0: minimize the norm of plain (u_k-u_{k-1})
%   u_ref: reference signal of u. It can be a constant vector or a function
%   handle that takes x as input.
%   DEBUG:
%   verbose_level: 
%       0: print no log (default)
%       1: print log of the rollout
%       2: print also the log of the controller
%   if 'verbose' is set to 1(0), it sets verbose_level = 1(0)
%% OUTPUTS
%   xs: trace of state (size: (plant_sys.xdim, total_k))
%   us: trace of control input (size: (control_sys.udim, total_k))
%       Note that the last control input in the array is actually not used
%       in the simulation. It is appended to match the size with xs and ts.
%   ts: trace of time (length: total_k)
%% Extra outputs (extraout)
%   Vs: trace of the CLF values (size: (control_sys.n_clf, total_k))
%   Bs: trace of the CLF values (size: (control_sys.n_cbf, total_k))
%   feas: trace of feasibility of the optimization-based controller (length: total_k)
%   slacks: trace of slack variables (size: (n_slack, total_k))
%       Note that n_slack depends on the choice of the controller.
%   comp_times: trace of computation time for each timestep (length: total_k)

settings = parse_function_args(varargin{:});
% Feedback-linearizer Determinant
if isa(control_sys, 'CtrlAffineSysFL')
    is_FL_system = true;
else
    is_FL_system = false;
end

if ~isfield(settings, 't0')
    t0 = 0;
else
    t0 = settings.t0;
end

% TODO: reverse calculate mu0 from u0
if ~isfield(settings, 'u0')
    u0 = [];
else
    if length(settings.u0) ~= control_sys.udim
        error("u0 dimension does not match with control_sys.udim.");
    end    
    u0 = settings.u0;
end

if ~isfield(settings, 'u_ref')
    u_ref = [];
else
    if isa(settings.u_ref, 'function_handle')
        error("Not supported in the current version.");
    elseif length(settings.u_ref) ~= control_sys.udim
        error("u_ref dimension does not match with control_sys.udim.");        
    end
    u_ref = settings.u_ref;
end

if ~isfield(settings, 'mu_ref')
    mu_ref = [];
else
    if isa(settings.mu_ref, 'function_handle')
        error("Not supported in the current version.");
    elseif length(settings.mu_ref) ~= control_sys.udim
        error("mu_ref dimension does not match with control_sys.udim.");        
    end
    mu_ref = settings.mu_ref;
end

if ~isfield(settings, 'dt')
    dt = 0.01;
else
    dt = settings.dt;
end

if ~isfield(settings, 'with_slack')
    with_slack = 1;
else
    with_slack = settings.with_slack;
end

if ~isfield(settings, 'ratio_u_diff')
    ratio_u_diff = 0;
else
    ratio_u_diff = settings.ratio_u_diff;
end

if ~isfield(settings, 'verbose_level')
    if isfield(settings, 'verbose')
        verbose_level = settings.verbose;
    else
        verbose_level = 0;
    end
else
    verbose_level = settings.verbose_level;
end

if ~isfield(settings, 'end_event')
    end_event = [];
else
    end_event = settings.end_event;
end

if ~isfield(settings, 'ode_func')
    ode_func = @ode45;
else
    ode_func = settings.ode_func;
end

% Dimension check
if plant_sys.udim ~= control_sys.udim
    error("plant_sys.udim and control_sys.udim should match each other.")
end
udim = plant_sys.udim;
if length(x0) ~= plant_sys.xdim
    error("x0 dimension does not match with plant_sys.xdim.");
end
if isrow(x0)
    x0 = x0';
end

% Initialize traces.
xs = x0;
ts = t0;
us = [];
% traces for extras, they remain empty if not returned.
feas = [];
Vs = [];
Bs = [];
feas = [];
slacks = [];
comp_times = [];

mus = [];
ys = [];
dys = [];

% Initialize state & time.
x = x0;
t = t0;

u_prev = zeros(control_sys.udim, 1);
mu_prev = zeros(control_sys.udim, 1);

end_simulation = false;
%% Run simulation.
% _t indicates variables for the current loop.
while ~end_simulation
    %% Determine control input.
    % TODO: can be simpler, but confusing.
    if t == t0 && ~isempty(u0)
        u = u0;
    else
        if isempty(u_ref)
            u_ref_t = ratio_u_diff * u_prev;
        else
            u_ref_t = u_ref;
        end
        if is_FL_system
            % TODO: going to support u_ref, mu_ref ver.
            %[u, extra_t]  = controller(x, u_ref_t, mu_ref_t, with_slack, verbose);
            % mu_ref_t should be supported
            if isempty(mu_ref)
                mu_ref_t = ratio_u_diff * u_prev;
            else
                mu_ref_t = mu_ref;
            end
            [u, extra_t] = controller(x, mu_ref_t, with_slack, verbose);
        else
            [u, extra_t] = controller(x, 'u_ref', u_ref_t, ...
                'with_slack', with_slack, 'verbose', (verbose_level>=2));
        end        
    end
    if verbose_level >= 1
        print_log(t, x, u, extra_t);
    end

    us = [us, u];
    feas = [feas, extra_t.feas];
    comp_times = [comp_times, comp_times];
    if with_slack
        slacks = [slacks, extra_t.slack];
    end
    if isfield(extra_t, 'Vs')
        Vs = [Vs, extra_t.Vs];
    end
    if isfield(extra_t, 'Bs')
        Bs = [Bs, extra_t.Bs];
    end
    
    %% Run simulation for one time step.
    t_end_t = min(t + dt, t0+T);
    if ~isempty(end_event)
        ode_opt = odeset('Events', end_event);
        [ts_t, xs_t, t_event] = ode_func( ...
            @(t, x) plant_sys.dynamics(t, x, u), ...
            [t, t_end_t], x, ode_opt);
        end_simulation = (ts_t(end) == t0 + T) || ~isempty(t_event);
    else
        [ts_t, xs_t] = ode_func(@(t, x) plant_sys.dynamics(t, x, u), ...
            [t, t_end_t], x);
        end_simulation = ts_t(end) == t0 + T;
    end            
    t = ts_t(end);
    x = xs_t(end, :)';
    u_prev = u;

    %% Record traces.
    xs = [xs, x];
    ts = [ts, t];
    if is_FL_system
        ys = [ys, extra_t.y];
        dys = [dys, extra_t.dy];
        mus = [mus, extra_t.mu];
    end
end % end of the main while loop
%% Add control input for the final timestep.
if isempty(u_ref)
    u_ref_t = ratio_u_diff * u_prev;
else
    u_ref_t = u_ref;
end
[u, extra_t] = controller(x, 'u_ref', u_ref_t, ...
                'with_slack', with_slack, 'verbose', (verbose_level>=2));
if verbose_level >= 1
    print_log(t, x, u, extra_t);
end

us = [us, u];
feas = [feas, extra_t.feas];
comp_times = [comp_times, comp_times];
if with_slack
    slacks = [slacks, extra_t.slack];
end
if isfield(extra_t, 'Vs')
    Vs = [Vs, extra_t.Vs];
end
if isfield(extra_t, 'Bs')
    Bs = [Bs, extra_t.Bs];
end
if is_FL_system
    ys = [ys, extra_t.y];
    dys = [dys, extra_t.dy];
    mus = [mus, extra_t.mu];
end

%% Make extraout
extraout.feas = feas;
extraout.comp_times = comp_times;
if ~isempty(slacks)
    extraout.slacks = slacks;
end
if ~isempty(Vs)
    extraout.Vs = Vs;
end
if ~isempty(Bs)
    extraout.Bs = Bs;
end
if is_FL_system
    extraout.mus = mus;
    extraout.ys = ys;
    extraout.dys = dys;
end
end % end of the main function.

function print_log(t, x, u, extra_t)
        fprintf("t: %.3f, \t x: ", t);
        fprintf("%.2g, ", x);
        fprintf("\t u: ");
        fprintf("%.2g, ", u);
        fprintf("\t feas: %d", extra_t.feas);
        % Add custom log here.
        fprintf("\n");
end
