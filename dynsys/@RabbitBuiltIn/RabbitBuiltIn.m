%% Main Reference
% Aaron Ames et al. Control Barrier Function based Quadratic Programs 
% with Application to Adaptive Cruise Control, CDC 2014, Table 1.

classdef RabbitBuiltIn < CtrlAffineSysFL
    properties
        reset_map_matrix = [1, 0, 0, 0, 0, 0, 0;
                             0, 1, 0, 0, 0, 0, 0;
                             0, 0, 1, 0, 0, 0, 0;
                             0, 0, 0, 0, 0, 1, 0;
                             0, 0, 0, 0, 0, 0, 1;
                             0, 0, 0, 1, 0, 0, 0;
                             0, 0, 0, 0, 1, 0, 0];
        fall_threshold = 0.699;
    end
    methods
        function obj = RabbitBuiltIn(params)
            % Always using built-in option for setup.
            obj = obj@CtrlAffineSysFL(params, 'built-in');
        end
        
        function [f, g] = defineDynamics(obj, x)
             n = obj.xdim / 2;

            q = x(1:n); 
            dq = x(n+1:2*n);

            %[D, C, G, B] = obj.gen_DCGB_rel(q);
            %D = D_gen_rel(q, obj.params.scale, obj.params.torso_add); % 7 x 7
            D = obj.D_gen_rel(q);
            C = obj.C_gen_rel(q,dq); % 7 x 7
            G = obj.G_gen_rel(q);
            B =[0 0 0 0 
                 0 0 0 0
                 0 0 0 0
                 1 0 0 0
                 0 1 0 0
                 0 0 1 0
                 0 0 0 1];
            JSt = J_RightToe(q);
            JSt = JSt([1,3],:);
            dJSt = dJ_RightToe(q,dq);
            dJSt = dJSt([1,3],:);

            H = C*dq + G;
            % Constraint Force to enforce the holonomic constraint:
            FSt_u = - pinv(JSt*(D\JSt'))*(JSt*(D\B));
            FSt_nu = - pinv(JSt*(D\JSt'))*(JSt*(D\(-H)) + dJSt*dq);

            f = [dq; D\(-C*dq - G)];
            g1 = [zeros(length(dq),obj.udim); D\B];
            g2 = [zeros(length(dq),size(FSt_nu ,1)); D\JSt'];

            % dx = f(x)+g(x)u of nominal model
            f = f + g2 * FSt_nu;
            g = g1 + g2 * FSt_u;
        end
        function f_ = f(obj, x)
            % RABBIT f(x)
            % TODO:multiple calculation of prior procedures => if we
            % maintain this architecture of different function, it is
            % inevitable or else, we can have another prior function but
            % still we need to waste memory space!
            % memory <-> time complexity issue
            [f_,~] = obj.defineDynamics(x);
        end
        function g_ = g(obj, x)
            % RABBIT g(x)
            [~, g_] = obj.defineDynamics(x);
        end
        function zs = z(obj, xs)
            % zero dynamics
            n = obj.xdim/2;
            q = xs(1:n, :);
            dq = xs(n+1:2*n, :);
            
            theta = q(3, :) + q(4, :) + q(5, :)/2;
%             theta_init = obj.params.legacy.theta_init;
%             theta_end = obj.params.legacy.theta_end;
%             phase = min(max((theta - theta_init)/(theta_end - theta_init)));
            
            dtheta = dq(3, :) + dq(4, :) + dq(5, :)/2;
            
            zs = [theta; dtheta];
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
            
            yd = obj.bezier_outputs(a_bez, phase); % Not very structured actually
            ya = q(4:end);
            y_ = ya-yd;
        end
        
%         function phase_ = phase(obj, x)
%             % TODO: Rabbit's phase => redundant(?)
%         end
        
        function Lfy = lf_y(obj, x)
            % TODO: Not used in RABBIT
            n = obj.xdim / 2;
            q = x(1:n);
            dq = x(n+1: 2*n);
            
            a_bez = obj.params.legacy.a_bez;
            theta_init = obj.params.legacy.theta_init;
            theta_end = obj.params.legacy.theta_end;
            
            dydq_des = obj.dydq_d_gen(q, a_bez, theta_init, theta_end); %[4 x 7]
            dydq_act = obj.dydq_a_gen(q);
            
            Lfy = (dydq_act - dydq_des)*dq; % Use dy = Lfy + Lgyu = Lfy property
        end
        
        function Lgy = lg_y(obj, x)
            % TODO: Not used in RABBIT
            Lgy = zeros(obj.ydim, 1);
        end
        
        function LgLfy = lglf_y(obj, x)
            % TODO: Rabbit's lie derivative of LgLfy
            n = obj.xdim / 2;
            q = x(1:n);
            
            a_bez = obj.params.legacy.a_bez;
            theta_init = obj.params.legacy.theta_init;
            theta_end = obj.params.legacy.theta_end;
            
            dydq_des = obj.dydq_d_gen(q, a_bez, theta_init, theta_end); %[4 x 7]
            dydq_act = obj.dydq_a_gen(q);
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
            
            L2ydq_des = obj.L2ydq_d_gen(q, dq, a_bez, theta_init, theta_end);
            L2ydq_act = obj.L2ydq_a_gen(q, dq);
            f_s = obj.f(x);
            
            L2fy = (L2ydq_act - L2ydq_des)*f_s;
        end
    end
end

