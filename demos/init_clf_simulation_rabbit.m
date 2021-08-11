%% Initialize Global Simulation Configuration of Rabbit
function dict = init_clf_simulation_rabbit
    %% Aggregate simulation setting in the params.config
    % TODO: we should classify some notions into certain groups, mixed ugly
    % now!
    
    %% Control Mode: You can use various control options
    dict.mode.motor_dynamics = 1; % going to be supported
    dict.mode.apply_rl_adaptation = 0; % going to be supported
    % 0: Original(Westervelt) 1: CLF-QP   2: L1-CLF-QP  3: RL-CLF-QP
    % IO-RL+CLF: apply_rl_adaptation=0, control_type=1
    dict.mode.control_type = 1;  % going to be supported
    
    %% Physical Property
    dict.physics.mu = 0.8; % friction constraint;
    
    %% Legacy : Trajectory optimization legacy, Initial setting
    % This is derived from the code Rabbit-IO-RL code
    % This is delicate, keep it fixed
    dict.legacy.theta_end = 3.384545423940834;
    dict.legacy.theta_init = 2.903323878938605;
    dict.legacy.a_bez = [2.42225522164138,2.50944978892455,2.66152682946601,2.78780950883388,2.88091497548752,2.94002391033518;0.667715484387163,0.638300798629787,0.490574595122909,0.478548440963253,0.473893838339675,0.576633115386526;2.94201463482336,2.99945422582175,2.85290119992043,2.52843337205465,2.42762552077273,2.44116379159792;0.575072495911993,0.635317622337650,1.13813693839107,1.22909640103678,1.00446628056820,0.717766128548961];
    dict.legacy.q0_ini = [-0.161463255651884;0.739782734886561;0.152681626367748;2.44507104015569;0.658106151373543;2.95111227794945;0.599243461678844];
    dict.legacy.dq0_ini = [0.781190436522269;0.243063912942872;0.247793330530586;1.10287289727466;-0.548512132037030;0.203061891957879;1.50818729243256];
    
    %% Fulfill params
    dict.xdim = length(dict.legacy.q0_ini) + length(dict.legacy.dq0_ini);
    dict.udim = 4; % TODO hard-coding
    dict.rel_deg_y = 2;
    dict.use_phase = false;
    dict.scale = 1;
    
    %% Input options
    dict.uSat = 150;
    dict.u_max = dict.uSat;
    dict.u_min  = -dict.uSat;

    %% Control Barrier Function Parameters
    dict.cbf.rate = 3;
    dict.cbf_gamma0 = 1;

    %% Control Lyapunov Function Parameters
    dict.clf.rate = 0.7;
    dict.clf.use_user_defined_rate = false; % set this to true if you don't want to use 
    
    dict.epsilon_FL = 0.045; % tune this for controlling convergence rate
    dict.Kp_FL = 1;
    dict.Kd_FL = 1.8;
    
    %% Quadratic Programming Parameters
    % objective function = A * |mu|^2 + B * |delta|^2
    % tune A and B to control the weight of each variable
    dict.weight.input = 1; % A (can be scalar or matrix)
    dict.weight.slack = 100; % B (should be scalar)
end