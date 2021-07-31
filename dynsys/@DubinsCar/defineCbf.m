function cbf = defineCbf(~, params, symbolic_state)
    x = symbolic_state;
    p_x = x(1); p_y = x(2); theta = x(3);

    v = params.v;
    xo = params.xo;
    yo = params.yo;
    d = params.d;
    if length(yo) ~= length(xo) || length(d) ~= length(xo)
        error("wrong obstacles' parameters")
    end
    
    n_obstacles = length(xo);
    cbf = cell(n_obstacles, 1);
    for i = 1:n_obstacles
        distance = (p_x - xo(i))^2 + (p_y - yo(i))^2 - d(i)^2;
        derivDistance = 2*(p_x-xo(i))*v*cos(theta) + 2*(p_y-yo(i))*v*sin(theta);
        cbf{i} = derivDistance + params.cbf_gamma0 * distance;
    end
end