classdef CtrlAffineSysFL < CtrlAffineSys
    %% Control-Affine Dynamic System with Feeback Linearization utilities
    properties
        % Autonomous matrix of the Linearized output dynamics
        F_FL
        F_FL_eps % auxilary
        % Actuation matrix of the Linearized output dynamics
        G_FL
        % Dimension of output
        ydim
        % Relative degree of y
        rel_deg_y
        % Output as a function handle (assuming relative degree 2)
        y
        phase
        % 1st order Lie derivative of the output
        lf_y
        lg_y
        % 2nd order Lie derivatives of the output as function handles
        lglf_y
        l2f_y
        % output related functions needed handle phase bounds
        y_max_exceed
        lf_y_max_exceed
        lg_y_max_exceed
        lglf_y_max_exceed
        l2f_y_max_exceed
        y_min_exceed
        lf_y_min_exceed
        lg_y_min_exceed
        lglf_y_min_exceed
        l2f_y_min_exceed        
        
        % CLF under feedback linearization and its derivatives as function handles.
        Gram_clf_FL % Gram matrix of the CLF for feedback linearization.
    end
    
    methods
        function obj = CtrlAffineSysFL(params)
            obj@CtrlAffineSys(params);
            
            [s, f, g] = defineSystem(obj, params);
            [y, phase, y_max_exceed, y_min_exceed] = obj.defineOutput(params, s);
            obj.initOutputDynamics(s, f, g, y, phase, y_max_exceed, y_min_exceed, params);
        end
        
        function [y, phase, y_max_exceed, y_min_exceed] = defineOutput(obj, params, symbolic_state)
            y = [];
            phase = [];
            y_max_exceed = [];
            y_min_exceed = [];
        end
        
        function V = clf_FL(obj, y, dy)
            eta_eps = [(1/obj.params.epsilon_FL)*y; dy];
            V = transpose(eta_eps)*obj.Gram_clf_FL*eta_eps;
        end
        
        function lF_clf_ = lF_clf_FL(obj, y, dy)
            eta_eps = [(1/obj.params.epsilon_FL)*y; dy];
            P = obj.Gram_clf_FL;
            lF_clf_ = transpose(eta_eps)*(transpose(obj.F_FL_eps)*P+P*obj.F_FL_eps)*eta_eps;
        end
        
        function lG_clf_ = lG_clf_FL(obj, y, dy)
            eta_eps = [(1/obj.params.epsilon_FL)*y; dy];
            P = obj.Gram_clf_FL;
            lG_clf_ = (2 * (obj.G_FL'*P) * eta_eps)';
        end
        
        function [y, dy, L2fy, LgLfy, phase] = eval_y(obj, s)
            if isempty(obj.phase)
                y = obj.y(s);
                dy = obj.lf_y(s);
                L2fy = obj.l2f_y(s);
                LgLfy = obj.lglf_y(s);
                phase = [];
                return
            end
            phase = obj.phase(s);
            if phase > obj.params.phase_max
                y = obj.y_max_exceed(s);
                dy = obj.lf_y_max_exceed(s);
                L2fy = obj.l2f_y_max_exceed(s);
                LgLfy = obj.lglf_y_max_exceed(s);
            elseif phase < obj.params.phase_min
                y = obj.y_min_exceed(s);
                dy = obj.lf_y_min_exceed(s);
                L2fy = obj.l2f_y_min_exceed(s);
                LgLfy = obj.lglf_y_min_exceed(s);
            else
                y = obj.y(s);
                dy = obj.lf_y(s);
                L2fy = obj.l2f_y(s);
                LgLfy = obj.lglf_y(s);
            end
        end
    end
end
