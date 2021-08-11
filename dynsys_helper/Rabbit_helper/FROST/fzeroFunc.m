function constr= fzeroFunc(z0, x0)
    
    yInit = z0(1);
    qSwingKneeInit = z0(2);
    dxInit = z0(3);
    dyInit = z0(4);
    
    qInit = x0(1:7);
    qInit(2) = yInit;
    qInit(7) = qSwingKneeInit;
    pFootStance = p_RightToe(qInit);
    pFootSwing = p_LeftToe(qInit);
    
    dqInit = x0(8:14);
    dqInit(1) = dxInit;
    dqInit(2) = dyInit;
    
    vFootStance = J_RightToe(qInit)*dqInit;
    vFootStance = vFootStance([1,3]);
    
    constr = [pFootStance(3);
              pFootSwing(3);
              vFootStance;
              ];
          

    
end