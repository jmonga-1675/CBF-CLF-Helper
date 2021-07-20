%% Author: Jason Choi (jason.choi@berkeley.edu)
function [mu, slack, B, V, feas] = ctrlCbfClfQpFL(obj, s, mu_ref, with_slack, verbose)
    %% Implementation of CBF-CLF-QP under feedback linearization structure.
    % Inputs:   s: state
    %           mu_ref: reference virtual control input
    %           with_slack: flag for relaxing the clf constraint(1: relax, 0: hard-constraint)
    %           verbose: flag for logging (1: print log, 0: run silently)
    % Outputs:  u: control input as a solution of the CBF-CLF-QP
    %           slack: slack variable for relaxation. (empty list when with_slack=0)
    %           B: Value of the CBF at current state.
    %           V: Value of the CLF at current state.
    %           feas: 1 if QP is feasible, 0 if infeasible. (Note: even
    %           when qp is infeasible, u is determined from quadprog.)
    if isempty(obj.cbf)
        error('CBF is not defined so ctrlCbfClfQpFL cannot be used. Create a class function [defineCbf] and set up cbf with symbolic expression.');
    end
        
    if nargin < 3
        mu_ref = zeros(obj.udim, 1);
    end
    if nargin < 4
        with_slack = 1;
    end
    if nargin < 5
        % Run QP without log in default condition.
        verbose = 0;
    end

    if size(mu_ref, 1) ~= obj.udim
        error("Wrong size of mu_ref, it should be (udim, 1) array.");
    end                
            
    V = obj.clf_FL(s);
    LfV = obj.lF_clf_FL(s);
    LgV = obj.lG_clf_FL(s);

    B = obj.cbf(s);
    LfB = obj.lf_cbf(s);
    LgB = obj.lg_cbf(s);
        
    y = obj.y(s);
    L2fy = obj.l2f_y(s);
    LgLfy = obj.lglf_y(s);
    inv_LgLfy = inv(LgLfy);
    u_star = -LgLfy\L2fy; %feedforward term
    if with_slack
        %% Decision variables (obj.udim + 2)
        %%      1-udim: \mu
        %%      udim+1: virtual cbf input (lie derivative)
        %%      udim+2: slack
        A = [LgV, 0, -1;
            zeros(1, obj.udim), -1, 0];
        b = [-LfV - obj.params.clf.rate * V;
            obj.params.cbf.rate * B];
        Aeq = [LgB * inv_LgLfy, -1, 0];
        beq = [-LfB-LgB*u_star];
        %% Add input constraints if u_max or u_min exists.
        if isfield(obj.params, 'u_max')
            A = [A; inv_LgLfy, zeros(obj.udim, 2);];
            if size(obj.params.u_max, 1) == 1
                b = [b; obj.params.u_max * ones(obj.udim, 1)-u_star];
            elseif size(obj.params.u_max, 1) == obj.udim
                b = [b; obj.params.u_max-ustar];
            else
                error("params.u_max should be either a scalar value or an (udim, 1) array.")
            end
        end
        if isfield(obj.params, 'u_min')
            A = [A; -inv_LgLfy, zeros(obj.udim, 2);];
            if size(obj.params.u_min, 1) == 1
                b = [b; u_star-obj.params.u_min * ones(obj.udim, 1)];
            elseif size(obj.params.u_min, 1) == obj.udim
                b = [b; u_star-obj.params.u_min];
            else
                error("params.u_min should be either a scalar value or an (udim, 1) array")
            end
        end        
    else
        %% Decision variables (obj.udim + 1)
        %%      1-udim: \mu
        %%      udim+1: virtual cbf input (lie derivative)
        A = [LgV, 0;
            zeros(1, obj.udim), -1];
        b = [-LfV - obj.params.clf.rate * V;
            obj.params.cbf.rate * B];
        Aeq = [LgB * inv_LgLfy, -1];
        beq = [-LfB-LgB*u_star];
        %% Add input constraints if u_max or u_min exists.
        if isfield(obj.params, 'u_max')
            A = [A; inv_LgLfy, zeros(obj.udim, 1);];
            if size(obj.params.u_max, 1) == 1
                b = [b; obj.params.u_max * ones(obj.udim, 1)-u_star];
            elseif size(obj.params.u_max, 1) == obj.udim
                b = [b; obj.params.u_max-ustar];
            else
                error("params.u_max should be either a scalar value or an (udim, 1) array.")
            end
        end
        if isfield(obj.params, 'u_min')
            A = [A; -inv_LgLfy, zeros(obj.udim, 1);];
            if size(obj.params.u_min, 1) == 1
                b = [b; u_star-obj.params.u_min * ones(obj.udim, 1)];
            elseif size(obj.params.u_min, 1) == obj.udim
                b = [b; u_star-obj.params.u_min];
            else
                error("params.u_min should be either a scalar value or an (udim, 1) array")
            end
        end
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
        options =  optimset('Display','notify');
    else
        options =  optimset('Display','off');
    end
    if with_slack         
        % cost = 0.5 [u' slack] H [u; slack] + f [u; slack]
        H = [weight_input, zeros(obj.udim, 1);
            zeros(1, obj.udim), obj.params.weight.slack];
        f_ = [-weight_input * mu_ref; 0];
        [mu_slack, ~, exitflag, ~] = quadprog(H, f_, A, b, Aeq, beq, [], [], [], options);
        if exitflag == -2
            feas = 0;
            disp("Infeasible QP. CBF constraint is conflicting with input constraints.");
        else
            feas = 1;
        end
        mu = mu_slack(1:obj.udim);
        slack = mu_slack(end);
    else
        % cost = 0.5 u' H u + f u
        H = weight_input;
        f_ = -weight_input * mu_ref;
        [mu, ~, exitflag, ~] = quadprog(H, f_, A, b, Aeq, beq, [], [], [], options);
        if exitflag == -2
            feas = 0;
            disp("Infeasible QP. CBF constraint is conflicting with input constraints.");
        else
            feas = 1;
        end
        slack = [];
    end
end