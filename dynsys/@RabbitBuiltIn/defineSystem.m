function [f, g] = defineSystem(obj, x)
% define f and g of the dynamical system
% Forked from RABBIT_IO_RL
% f and g would be calculated - duplicated
% => TODO: should be incorporated
	n = obj.xdim / 2;
    if obj.params.mode.motor_dynamics
        u_s = x(end-3 : end);
    end
    
    q = x(1:n); 
    dq = x(n+1:2*n);
    
    Dm = D_gen_rel(q,1); % 7 x 7
    Cm = C_gen_rel(q,dq,1); % 7 x 7
    Gm = G_gen_rel(q,1);
    Bm =[0 0 0 0 
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
    
    Hm= Cm*dq + Gm;
    % Constraint Force to enforce the holonomic constraint:
    FSt_u_m = - pinv(JSt*(Dm\JSt'))*(JSt*(Dm\Bm));
    FSt_nu_m = - pinv(JSt*(Dm\JSt'))*(JSt*(Dm\(-Hm)) + dJSt*dq);

    f_m = [dq; Dm\(-Cm*dq - Gm)];
    g1_m = [zeros(length(dq),obj.udim); Dm\Bm];
    g2_m = [zeros(length(dq),size(FSt_nu_m,1)); Dm\JSt'];

    % dx = f(x)+g(x)u of nominal model
    f = f_m + g2_m*(FSt_nu_m);
    g = g1_m + g2_m*FSt_u_m;
end

