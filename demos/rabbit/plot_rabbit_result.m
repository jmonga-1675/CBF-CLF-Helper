function plot_rabbit_result(xs, ts, us, extras)
% Plot the result of the RABBIT model
% Argument
%   xs: trajectory of RABBIT simulation
%   ts: time stamp of RABBIT simulation
%   us: control input of RABBIT simulation
%   extras: extra logs of RABBIT simulation

% Unravel extras data into one array
Vs = []; dVs_error = []; etas_zs = [];
feas = []; slacks = []; mus = []; ys = []; dys = []; Fsts = [];
for i = 1:length(extras)
    Vs = [Vs; extras(i).logs.Vs];
    dVs_error = [dVs_error; extras(i).logs.dVs_error'];
    etas_zs = [etas_zs; extras(i).logs.etas_zs];
    dys = [dys; extras(i).dys];
    ys = [ys; extras(i).ys];
    
    feas = [feas; extras(i).feas];
    slacks = [slacks; extras(i).slacks];
    mus = [mus; extras(i).mus];
    Fsts = [Fsts; extras(i).Fsts];
end

% Plot setting
subplot_gap = [0 0.1];
fig_sz = [8 3]; 
fig_sz_with_subplots = [8 4]; 
plot_pos = [0 0 8 3];

main_color_orange = [217 83 25]/255; % r,g,b
main_color_dark_red = [162 20 47]/255; % r,g,b
main_color_blue = [0 114 189]/255; % r,g,b
main_color_green = [119 172 48]/255; 

figure(1);
plot(ts, ys(:,1)); grid on;
xlabel('Time(s)'); ylabel('q_1');
title('Stance Leg Hip Angle');

figure(2);
plot(ts, ys(:,2)); grid on;
xlabel('Time(s)'); ylabel('q_3');
title('Stance Leg Knee Angle');

figure(3);
plot(ts, ys(:,3)); grid on;
xlabel('Time(s)'); ylabel('q_2');
title('Swing Leg Hip Angle');

figure(4);
plot(ts, ys(:,4)); grid on;
xlabel('Time(s)'); ylabel('q_4');
title('Swing Leg Knee Angle');

figure(5); plot(ts, Fsts(:,1)); grid on; hold on; xlabel('t');
ylabel('FSt'); plot(ts, Fsts(:,2)); legend(['FSt(1)';'FSt(2)']);

figure(6);
plot(ts,abs(Fsts(:,1)./Fsts(:,2)), 'color', main_color_blue, 'LineWidth', 1); grid on; hold on;
title('Contact Force Ratio ($|F_{T}$/$F_{N}|$)', 'Interpreter', 'latex');
xlim([0 ts(end)+1]);
xlabel('Time (sec)');
ylabel('$|F_{T}$/$F_{N}|$', 'Interpreter', 'latex');
callPlotSettings(fig_sz, plot_pos);
%print('-f6','plot_F_ratio', '-dpng');

figure(7);
plot(ts, us); grid on; hold on;
title('torques'); legend('u1','u2','u3','u4');

figure(9); 
subplot1(2,1,'Gap', subplot_gap);
subplot1(1); hold on;
title('Euclidian Norm of the Tracking Error');
plot(ts,sqrt(ys(:,1).^2+ys(:,2).^2+...
    ys(:,3).^2+ys(:,4).^2), 'color', main_color_blue); 
ylabel('||y||_2 (rad)'); grid on;
callPlotSettings(fig_sz_with_subplots, plot_pos);

subplot1(2); hold on;
title('Euclidian Norm of the Derivative of the Tracking Error');
plot(ts,sqrt(dys(:,1).^2+dys(:,2).^2+...
    dys(:,3).^2+dys(:,4).^2), 'color', main_color_blue); 
ylabel('$||\dot{y}||_{2}$ (rad/s)','Interpreter','latex'); xlabel('Time (s)'); grid on;
callPlotSettings(fig_sz_with_subplots, plot_pos);
%print('-f9','plot_tracking_error', '-dpng');
end

