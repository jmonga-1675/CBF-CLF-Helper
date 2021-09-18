function [u, y, feedforward, u_raw] = ctrlFeedbackLinearize(obj, s, mu)
%% Implementation of feedback linearization structure
% Inputs:   s: state
%           mu: virtual control input
% Outputs:  u: control input.
% u = u_star + LgLf\mu
y = obj.y(s);
if obj.rel_deg_y == 1
    Lfy = obj.lf_y(s);
    Lgy = obj.lg_y(s);
    if rank(Lgy) < obj.ydim
        error("Feedback Linearization fail because lg_y is not invertible")
    end
    feedforward = -Lgy\Lfy;
    u_raw = feedforward + Lgy\mu;
    
elseif obj.rel_deg_y == 2
    LgLfy = obj.lglf_y(s);
    L2f_y = obj.l2f_y(s);
    if rank(LgLfy) < obj.ydim
        error("Feedback Linearization fail because lglf_y is not invertible")
    end
    feedforward = -LgLfy\L2f_y;
    u_raw = feedforward + LgLfy\mu;    
end

u = u_raw;

%% Clip input into the physical constraint
eps = 1e-4;
if isfield(obj.params, 'u_max')
    if size(obj.params.u_max, 1) == 1 % row vector
        for i = 1:obj.udim
            if u(i) > obj.params.u_max + eps
                u(i) = obj.params.u_max;
                disp("Warning: virtual control input does not satisfy input constraint.")
            end
        end
    elseif size(obj.params.u_max, 1) == obj.udim % column vector
        for i = 1:obj.udim
            if u(i) > obj.params.u_max(i) + eps
                u(i) = obj.params.u_max(i);
                disp("Warning: virtual control input does not satisfy input constraint.")
            end    
        end
    else
        error("params.u_max Unknown size.")
    end    
end

if isfield(obj.params, 'u_min')
    if size(obj.params.u_min, 1) == 1
        for i = 1:obj.udim
            if u(i) < obj.params.u_min - eps
                u(i) = obj.params.u_min;
                disp("Warning: virtual control input does not satisfy input constraint.")
            end
        end
    elseif size(obj.params.u_min, 1) == obj.udim
        for i = 1:obj.udim
            if u(i) < obj.params.u_min(i) - eps
                u(i) = obj.params.u_min(i);
                disp("Warning: virtual control input does not satisfy input constraint.")
            end    
        end
    else
        error("params.u_min Unknown size.")
    end    
end