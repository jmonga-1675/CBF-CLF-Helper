function clf = defineClf(obj, params, symbolic_state)
    s = symbolic_state;
    E = getEnergy(obj, s);
    E_ref = params.m * params.g * params.l;

    l = obj.params.l;  % [m]      length of pendulum
    m = obj.params.m;  % [kg]     mass of pendulum
    M = obj.params.M;  % [kg]     mass of cart
    b = obj.params.b;  % [N/m/s]  coefficient of friction between cart and ground
    g = obj.params.g; % [m/s^2]  acceleration of gravity
    A = [0, 1, 0, 0;
        0, 0, -m*g/M, 0;
        0, 0, 0, 1;
        0, 0, (m+M)*g/M, 0];
    B = [0; 1/M; 0; -1/(M*l)];
    Q = eye(4);
    Q(1, 1) = 0.1;
    Q(2, 2) = 0.1;
    Q(3, 3) = 100;
    Q(4, 4) = 20;
    R = 0.1;
    [K, P] = lqr(A, B, Q, R);
    disp(P)
    disp(eig(A-B*K))
    s_n = s;
    s_n(3) = sin(s(3));
    s_n(4) = -s(4);
    V_neighbor = s_n' * P * s_n;

    % clf =  (1-cos(s(3))) * V_neighbor;
    % clf =  V_neighbor;
    % clf = (E - E_ref)^2;
    % clf = (E - E_ref)^2 + 5 * V_neighbor;
    clf = (E - E_ref)^2 + 30 * (1-cos(s(3))) * V_neighbor;
    % clf = abs(E - E_ref) * V_neighbor;

    % clf = (E - E_ref)^2 + 20 * exp(-(E - E_ref)^2) * V_neighbor;

    % clf = V_neighbor;
end