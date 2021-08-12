function [FSt_u, FSt_nu] = get_force(obj, x)
    n = obj.xdim / 2;
    if obj.params.mode.motor_dynamics
        u_s = x(end-3 : end);
    end
    
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
end

