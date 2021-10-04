%% Author: Jason Choi (jason.choi@berkeley.edu)
function [u, extraout] = ctrlClfQp(obj, x, u_ref, with_slack, verbose)
    %% Implementation of the vanilla CLF-QP
    % Inputs:   x: state
    %           u_ref: reference control input
    %           with_slack: flag for relaxing (1: relax, 0: hard CLF constraint)
    %           verbose: flag for logging (1: print log, 0: run silently)
    % Outputs:  u: control input as a solution of the CLF-QP
    %   extraout:
    %           slack: slack variable for relaxation. (empty when with_slack=0)
    %           Vs: CLF values at the current state.
    %           feas: 1 if QP is feasible, 0 if infeasible. (Note: even
    %           when qp is infeasible, u is determined from quadprog.)
    %           comp_time: computation time to run the solver.

    if obj.n_clf == 0
        error('CLF is not set up so ctrlClfQp cannot be used.');
    end

    if nargin < 3 || isempty(u_ref)
        % If u_ref is given, CLF-QP minimizes the norm of u-u_ref        
        % Default reference control input is u.
        u_ref = zeros(obj.udim, 1);
    end
    if nargin < 4
        % Relaxing is activated in default condition.
        with_slack = 1;
    end
    if nargin < 5
        % Run QP without log in default condition.
        verbose = 0;
    end
    
    if size(u_ref, 1) ~= obj.udim
        error("Wrong size of u_ref, it should be (udim, 1) array.");
    end
    
    tstart = tic;
    Vs = obj.clf(x);
    % Lie derivatives of the CLF.
    LfVs = obj.lf_clf(x);
    LgVs = obj.lg_clf(x);

    %% Constraints : A[u; slack] <= b
    A = LgVs;
    b = -LfVs - obj.clf_rate * Vs;
    if ~isempty(obj.u_max)
        A = [A; eye(obj.udim)];
        b = [b; obj.u_max];
    end
    if ~isempty(obj.u_min)
        A = [A; -eye(obj.udim)];
        b = [b; -obj.u_min];
    end
    if with_slack
        A_slack = -eye(obj.n_clf);
        if ~isempty(obj.u_max)
            A_slack = [A_slack; zeros(obj.udim, obj.n_clf)];
        end
        if ~isempty(obj.u_min)
            A_slack = [A_slack; zeros(obj.udim, obj.n_clf)];
        end
        A = [A, A_slack];
    end
    
    %% Cost
    if isfield(obj.params.weight, 'input')
        if size(obj.params.weight.input, 1) == 1 
            weight_input = obj.params.weight.input * eye(obj.udim);
        elseif all(size(obj.params.weight.input) == obj.udim)
            weight_input = obj.params.weight.input;
        else
            error("params.weight.input should be either a scalar value or an (udim, udim) array.")
        end
    else
        weight_input = eye(obj.udim);
    end
    
    if verbose
        options =  optimoptions('quadprog', 'ConstraintTolerance', 1e-6, 'StepTolerance', 1e-12, 'Display','iter');
    else
        options =  optimoptions('quadprog', 'ConstraintTolerance', 1e-6, 'Display','off');
    end
    
    if with_slack         
        % cost = 0.5 [u' slack] H [u; slack] + f [u; slack]
        H = [weight_input, zeros(obj.udim, obj.n_clf);
            zeros(obj.n_clf, obj.udim), diag(obj.params.weight.slack)];
        f_ = [-weight_input * u_ref; zeros(obj.n_clf, 1)];
        [u_slack, ~, exitflag, ~] = quadprog(H, f_, A, b, [], [], [], [], [], options);
        if exitflag == -2            
            feas = 0;
            disp("Infeasible QP. Numerical error might have occured.");
            % Making up best-effort heuristic solution.            
            u = zeros(obj.udim, 1);
            for i = 1:obj.udim
                u(i) = obj.u_min(i) * (LgVs(i) > 0) + obj.u_max(i) * (LgVs(i) <= 0);
            end
            slack = zeros(obj.n_clf, 1);
        else
            feas = 1;
            u = u_slack(1:obj.udim);
            slack = u_slack(obj.udim+1:end);
        end
    else
        H = weight_input;
        f_ = -weight_input * u_ref;
        [u, ~, exitflag, ~] = quadprog(H, f_, A, b, [], [], [], [], [], options);
        if exitflag == -2
            feas = 0;
            disp("Infeasible QP. CLF constraint is conflicting with input constraints.");
            % Making up best-effort heuristic solution.
            u = zeros(obj.udim, 1);
            for i = 1:obj.udim
                u(i) = obj.u_min(i) * (LgVs(i) > 0) + obj.u_max(i) * (LgVs(i) <= 0);
            end
        else
            feas = 1;
        end
        slack = [];
    end
    comp_time = toc(tstart);
    
    extraout.slack = slack;
    extraout.Vs = Vs;
    extraout.feas = feas;
    extraout.comp_time = comp_time;
end
