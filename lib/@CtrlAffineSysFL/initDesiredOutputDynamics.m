function initDesiredOutputDynamics(obj, params)
%% set up desired linear output dynamics
    if ~isfield(params, 'epsilon_FL')
        error("Define variable in the params called 'epsilon_FL', which is the rate of RES-CLF");
    end
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
    eps = params.epsilon_FL;
    Kp = params.Kp_FL;
    Kd = params.Kd_FL;
    A = [zeros(obj.ydim), eye(obj.ydim);
        -Kp*eye(obj.ydim), -Kd*eye(obj.ydim)];
    P = lyap(A', Q);
    
    if isfield(params, 'clf')
        if ~isfield(params.clf, 'rate')
            obj.params.clf.rate = min(eig(Q))/max(eig(P));
            fprintf("Setting clf rate automatically to %.3f \n", obj.params.clf.rate);
        end
    else
        obj.params.clf.rate = min(eig(Q))/max(eig(P));
        fprintf("Setting clf rate automatically to %.3f \n", obj.params.clf.rate);
    end

    obj.F_FL = [zeros(obj.ydim), 1/eps * eye(obj.ydim);
        zeros(obj.ydim), zeros(obj.ydim)];
    obj.G_FL = [zeros(obj.ydim);
                eye(obj.ydim)];
    obj.clf = @(s) control_lyapunov_function(s, obj, eps, P);
    obj.lf_clf = @(s) lf_clf(s, obj, eps, P);
    obj.lg_clf = @(s) lg_clf(s, obj, eps, P);
end

function V = control_lyapunov_function(s, obj, eps, P)
    y = obj.y(s);
    dy = obj.lf_y(s);
    eta_eps = [(1/eps)*y; dy];
    V = transpose(eta_eps)*P*eta_eps;
end

function lf_clf_ = lf_clf(s, obj, eps, P)
    y = obj.y(s);
    dy = obj.lf_y(s);
    eta_eps = [(1/eps)*y; dy];
    lf_clf_ = transpose(eta_eps)*(transpose(obj.F_FL)*P+P*obj.F_FL)*eta_eps;
end

function lg_clf_ = lg_clf(s, obj, eps, P)
    y = obj.y(s);
    dy = obj.lf_y(s);
    eta_eps = [(1/eps)*y; dy];
    lg_clf_ = 2 * (obj.G_FL'*P) * eta_eps;
end

