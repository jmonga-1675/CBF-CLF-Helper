classdef CartPole < CtrlAffineSysFL
    methods
        function E = getEnergy(obj, s)
            l = obj.params.l;  % [m]      length of pendulum
            m = obj.params.m;  % [kg]     mass of pendulum
            M = obj.params.M;  % [kg]     mass of cart
            b = obj.params.b;  % [N/m/s]  coefficient of friction between cart and ground
            g = obj.params.g; % [m/s^2]  acceleration of gravity
            
            % kinetic energy
            T = 0.5 * (M+m) * s(2)^2 + m * s(2) * s(4) * l * cos(s(3)) + 0.5 * m * l^2 * s(4)^2;
            % potential energy
            U = -m * g * l * cos(s(3));
            E = T + U;
        end
    end
end
