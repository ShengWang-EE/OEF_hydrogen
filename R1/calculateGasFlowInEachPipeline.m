function gasFlow = calculateGasFlowInEachPipeline(auxiliaryVar, S_sqrt,composition_hy,composition_gas,signGf, mpc,gtd,nGl,M_hy,M_gas,M_air,R_air,T_stp,Prs_stp,T_gas,Z)
% auxiliaryVar = auxiliaryVar * 1e10; %Pa
for  m = 1:nGl
    if signGf(m) == 1 % positive
        fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);
    else
        fb = mpc.Gline(m,2); tb = mpc.Gline(m,1);
    end
    % para
    D = gtd.Gline(m,3); F = gtd.Gline(m,5);
    L = gtd.Gline(m,4); 

    %
%     gasFlow(m) = auxiliaryVar(m)/S_sqrt(fb) * ( (pi^2*R_air/64) * (T_stp/Prs_stp)^2 * ...
%         ( D^5) / (F*L*T_gas*Z) ) ^0.5;
    S_gas = M_gas/M_air;
    C = mpc.Gline(m,3);
    newC = C*sqrt(S_gas)*sqrt(Z_gas)/ S_sqrt(fb) / ;
    gasFlow(m) = newC*auxiliaryVar(m) * 1e6/24/3600; %m3/s
end
end