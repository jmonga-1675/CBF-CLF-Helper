function [s, f, g] = defineSystem(obj, params)

l = params.l;  % [m]      length of pendulum
m = params.m;  % [kg]     mass of pendulum
M = params.M;  % [kg]     mass of cart
b = params.b;  % [N/m/s]  coefficient of friction between cart and ground
g = params.g; % [m/s^2]  acceleration of gravity

syms x dx theta dtheta
% Theta: upright is pi (-pi)
s = [x; dx; theta; dtheta];
    
f = [s(2);
    ( 2*m*l*s(4)^2*sin(s(3)) + 3*m*g*sin(s(3))*cos(s(3)) ...
      - 4*b*s(2) )/( 4*(M+m)-3*m*cos(s(3))^2 );
    s(4);
    (-3*m*l*s(4)^2*sin(s(3))*cos(s(3)) - 6*(M+m)*g*sin(s(3)) ...
      - 6*(-b*s(2))*cos(s(3)) )/( 4*l*(m+M)-3*m*l*cos(s(3))^2 )];

g = [0;
    4/( 4*(M+m)-3*m*cos(s(3))^2 );
    0;
    (- 6*cos(s(3)) )/( 4*l*(m+M)-3*m*l*cos(s(3))^2 )];
