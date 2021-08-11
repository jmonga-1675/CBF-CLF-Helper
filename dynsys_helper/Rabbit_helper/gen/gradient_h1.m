function gradient_h1 = gradient_h1(in1,L3,L4,R1,l_max)
%GRADIENT_H1
%    GRADIENT_H1 = GRADIENT_H1(IN1,L3,L4,R1,L_MAX)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    24-Jan-2020 08:14:27

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = cos(q4);
t6 = sin(q1);
t7 = sin(q2);
t8 = sin(q3);
t9 = sin(q4);
t10 = L3.*t2;
t11 = L3.*t3;
t12 = L4.*t4;
t13 = L4.*t5;
t14 = L3.*t6;
t15 = L3.*t7;
t16 = L4.*t8;
t17 = L4.*t9;
t18 = -t11;
t19 = -t13;
t20 = -t15;
t21 = -t17;
t22 = t10+t12+t18+t19;
t23 = R1+t14+t16+t20+t21;
t24 = t22.^2;
t25 = t23.^2;
t26 = t24+t25;
t27 = 1.0./sqrt(t26);
gradient_h1 = [t27.*(t10.*t23.*2.0-t14.*t22.*2.0).*(-1.0./2.0);(t27.*(t11.*t23.*2.0-t15.*t22.*2.0))./2.0;t27.*(t12.*t23.*2.0-t16.*t22.*2.0).*(-1.0./2.0);(t27.*(t13.*t23.*2.0-t17.*t22.*2.0))./2.0;0.0];
