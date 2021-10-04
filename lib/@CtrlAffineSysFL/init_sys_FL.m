function init_sys_FL(obj, params)
%% Functions that initialize dynamic system
    if ~isfield(params, 'use_phase')
        obj.use_phase = 1;    
    else
        obj.use_phase = params.use_phase;
    end

    %% Symbolic Function
    if strcmp(obj.setup_option, 'symbolic')
        disp(['Setting up feedback linearization dynamics, CLFs, CBFs from defined symbolic expressions.', ...
            '(This might take time.)']);
        [x, f, g] = obj.defineSystem(params);
        [y, phase, y_max_exceed, y_min_exceed] = obj.defineOutput(params, x);
        obj.initOutputDynamics(x, f, g, y, phase, y_max_exceed, y_min_exceed, params);
    
    %% Non Symbolic Function
    elseif strcmp(obj.setup_option, 'built-in')
        if ~isfield(params, 'rel_deg_y')
            error("rel_deg_y should be specified for built-in setup.");
        end
        obj.rel_deg_y = params.rel_deg_y; % TODO: other way to judge it?
        obj.ydim = obj.udim;
    else
        error("Undefined setup_option.");
    end
    %% Set up desired linear output dynamics
    %% Both Symbolic and Non-symbolic
    if ~isfield(params, 'epsilon_FL')
        error("Define variable in the params called 'epsilon_FL', which is the rate of RES-CLF");
    end
    eps = params.epsilon_FL;
    obj.F_FL = [zeros(obj.ydim), eye(obj.ydim);
        zeros(obj.ydim), zeros(obj.ydim)];
    obj.F_FL_eps = [zeros(obj.ydim), 1/eps * eye(obj.ydim);
        zeros(obj.ydim), zeros(obj.ydim)];
    obj.G_FL = [zeros(obj.ydim);
                eye(obj.ydim)];

    %% Set up CLF.
    % In controlSysFL, CLF is automatically set, so obj.n_clf can be
    % hard-coded to 1
    obj.n_clf = 1;
    if ~isfield(params, 'Kp_FL')
        error("Define variable in the params called 'Kp_FL', p gain for RES-CLF");
    end
    if ~isfield(params, 'Kd_FL')
        error("Define variable in the params called 'Kd_FL', d gain for RES-CLF");
    end
    if ~isfield(params, 'Q_FL')
        disp("Setting Q matrix for CLF to default identity matrix.");
        Q = eye(2*obj.ydim);
    else
        if size(params.Q_FL, 1) ~= 2*obj.ydim
            error("Wrong Q_FL size.")
        end
        Q = params.Q_FL;
    end
    Kp = params.Kp_FL;
    Kd = params.Kd_FL;
    A = [zeros(obj.ydim), eye(obj.ydim);
        -Kp*eye(obj.ydim), -Kd*eye(obj.ydim)];
    obj.Gram_clf_FL = lyap(A', Q); % A'P+PA=Q, find P 

    if isfield(params, 'clf')
        if ~isfield(params.clf, 'rate') || ~params.clf.use_user_defined_rate
            obj.params.clf.rate = min(eig(Q))/max(eig(obj.Gram_clf_FL))/eps;
            fprintf("Setting clf rate automatically to %.3f \n", obj.params.clf.rate);
        end
    else
        obj.params.clf.rate = min(eig(Q))/max(eig(obj.Gram_clf_FL))/eps;
        fprintf("Setting clf rate automatically to %.3f \n", obj.params.clf.rate);
    end
end

