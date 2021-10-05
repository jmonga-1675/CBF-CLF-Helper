function x_after = reset_map(obj, x_before)
    n = obj.xdim/2;
    q = obj.reset_map_matrix * x_before(1:n)';
    dq = obj.reset_map_matrix * obj.dq_pos_gen_v2([x_before(1:n) x_before(n+1:2*n)]);
    x_after = [q; dq];
end