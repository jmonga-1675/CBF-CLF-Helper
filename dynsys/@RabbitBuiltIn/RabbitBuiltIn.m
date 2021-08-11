%% Main Reference
% Aaron Ames et al. Control Barrier Function based Quadratic Programs 
% with Application to Adaptive Cruise Control, CDC 2014, Table 1.

classdef RabbitBuiltIn < CtrlAffineSysFL    
    methods
        function obj = RabbitBuiltIn(params)
            % Always using built-in option for setup.
            obj = obj@CtrlAffineSysFL(params, 'built-in');
        end
        function f_ = f(obj, x)
            % RABBIT f(x)
            % TODO:multiple calculation of prior procedures => if we
            % maintain this architecture of different function, it is
            % inevitable or else, we can have another prior function but
            % still we need to waste memory space!
            % memory <-> time complexity issue
            [f_,~] = obj.defineSystem(x);
        end
        function g_ = g(obj, x)
            % RABBIT g(x)
            [~, g_] = obj.defineSystem(x);
        end
        function V = clf(obj, x)
            % TODO: Lyapunov Function: redundant actually
            n = obj.xdim/2;
            y = obj.y(x);
            dy = obj.lf_y(x);
            eta_eps = [(1/obj.params.epsilon_FL)*y; dy];
            P = obj.Gram_clf_FL;
            % In case of sanity check 
            if isempty(obj.Gram_clf_FL)
                P = eye(obj.udim * 2); % exploiting obj.udim = obj.ydim
            end
            V = transpose(eta_eps)*P*eta_eps;
        end
        function LfV = lf_clf(obj, x)
            % TODO: LfV(x) : redundant => CtrlAffineSysFL rehandled
            n = obj.xdim/2;
            y = obj.y(x);
            dy = obj.lf_y(x);
            eta_eps = [(1/obj.params.epsilon_FL)*y; dy];
            
            P = obj.Gram_clf_FL;
            F_FL_eps = obj.F_FL_eps;
            
            % In case of sanity check
            if isempty(obj.F_FL_eps)
                F_FL_eps = eye(2 * obj.udim);
                P = eye(2 * obj.udim);
            end
            LfV = transpose(eta_eps)*(transpose(F_FL_eps)*P+P*F_FL_eps)*eta_eps;
        end
        function LgV = lg_clf(obj, x)
            % TODO: LgV(x) : redundant => CtrlAffineSysFL rehandled
            n = obj.xdim/2;
            y = obj.y(x);
            dy = obj.lf_y(x);
            eta_eps = [(1/obj.params.epsilon_FL)*y; dy];
            
            P = obj.Gram_clf_FL;
            G_FL = obj.G_FL;
            
            if isempty(obj.G_FL) % In case of sanity check
                G_FL = [zeros(obj.udim);
                    eye(obj.udim)];
                P = eye(2*obj.udim);
            end
            
            LgV = (2 * (G_FL'*P) * eta_eps)';
        end
        function B = cbf(obj, x)
            % TODO: CBF(x)
        end
        function LfB = lf_cbf(obj, x)
            % TODO: LfB(x)
        end
        function LgB = lg_cbf(obj, x)
            % TODO: LgB(x)
        end
        
        function y_ = y(obj, x)
            % TODO: Rabbit's output
            n = obj.xdim / 2;
            q = x(1:n);
            theta = q(3) + q(4) + q(5)/2;
            
            a_bez = obj.params.legacy.a_bez;
            theta_init = obj.params.legacy.theta_init;
            theta_end = obj.params.legacy.theta_end;
            phase = min(max((theta - theta_init)/(theta_end - theta_init)));
            
            yd = bezier_outputs(a_bez, phase); % Not very structured actually
            ya = q(4:end);
            y_ = ya-yd;
        end
        
        function phase_ = phase(obj, x)
            % TODO: Rabbit's phase => redundant(?)
        end
        
        function Lfy = lf_y(obj, x)
            % TODO: Not used in RABBIT
            n = obj.xdim / 2;
            q = x(1:n);
            dq = x(n+1: 2*n);
            
            a_bez = obj.params.legacy.a_bez;
            theta_init = obj.params.legacy.theta_init;
            theta_end = obj.params.legacy.theta_end;
            
            dydq_des = dydq_d_gen(q, a_bez, theta_init, theta_end); %[4 x 7]
            dydq_act = dydq_a_gen(q);
            
            Lfy = (dydq_act - dydq_des)*dq; % Use dy = Lfy + Lgyu = Lfy property
        end
        
        function Lgy = lg_y(obj, x)
            % TODO: Not used in RABBIT
            Lgy = -1;
        end
        
        function LgLfy = lglf_y(obj, x)
            % TODO: Rabbit's lie derivative of LgLfy
            n = obj.xdim / 2;
            q = x(1:n);
            
            a_bez = obj.params.legacy.a_bez;
            theta_init = obj.params.legacy.theta_init;
            theta_end = obj.params.legacy.theta_end;
            
            dydq_des = dydq_d_gen(q, a_bez, theta_init, theta_end); %[4 x 7]
            dydq_act = dydq_a_gen(q);
            g_s = obj.g(x);
            
            LgLfy = (dydq_act-dydq_des)*g_s(n+1:end, :);
        end
        
        function L2fy = l2f_y(obj, x)
            % TODO: Rabbit's lie derivative of L2fy
            n = obj.xdim/2;
            q = x(1:n);
            dq = x(n+1:2*n);
            
            a_bez = obj.params.legacy.a_bez;
            theta_init = obj.params.legacy.theta_init;
            theta_end = obj.params.legacy.theta_end;
            
            L2ydq_des = L2ydq_d_gen(q, dq, a_bez, theta_init, theta_end);
            L2ydq_act = L2ydq_a_gen(q, dq);
            f_s = obj.f(x);
            
            L2fy = (L2ydq_act - L2ydq_des)*f_s;
        end
        
        function y_max_exceed_ = y_max_exceed(obj, x)
        end
        
        function lf_y_max_exceed_ = lf_y_max_exceed(obj, x)
        end
        
        function lg_y_max_exceed_ = lg_y_max_exceed(obj, x)
        end
        
        function lglf_y_max_exceed_ = lglf_y_max_exceed(obj, x)
        end
        
        function l2f_y_max_exceed_ = l2f_y_max_exceed(obj, x)
        end
        
        function y_min_exceed_ = y_min_exceed(obj, x)
        end
        
        function lf_y_min_exceed_ = lf_y_min_exceed(obj, x)
        end
        
        function lg_y_min_exceed_ = lg_y_min_exceed(obj, x)
        end
        
        function lglf_y_min_exceed_ = lglf_y_min_exceed(obj, x)
        end
        
        function l2f_y_min_exceed_ = l2f_y_min_exceed(obj, x)
        end
    end
end

