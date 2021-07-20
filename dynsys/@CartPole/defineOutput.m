function [y, phase, y_max_exceed, y_min_exceed] = defineOutput(obj, params, symbolic_state)
    x = symbolic_state(1);
    theta = symbolic_state(3);
%    y = params.l * sin(theta);    
    phase = x + 0.5 * theta;
    c_normal = 2.37079632679490; % normalization constant.
    phase = phase / c_normal;
    
    y_desired = obj.bezier(params.bezier, phase);
    y_desired_max = obj.bezier(params.bezier, params.phase_max);
    y_desired_min = obj.bezier(params.bezier, params.phase_min);
%     M = length(params.bezier)-1; % Degree of the bezier curve
%     y_desired = sym(0);
%     for i = 0:M
%         y_desired = y_desired + ...
%             params.bezier(i+1) * nchoosek(M, i) * (1-phase_var)^(M-i) * phase_var^i;
%     end
    
    y_actual = theta;
    y = y_actual - y_desired;
    y_max_exceed = y_actual - y_desired_max;
    y_min_exceed = y_actual - y_desired_min;    
end