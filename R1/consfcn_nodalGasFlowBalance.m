function f = consfcn_nodalGasFlowBalance(PGs,Qd,Qgpp,Qptg, gasFlow,mpc,nGasType,nGPP,nGasLoad,iGasLoad)
%% parameter
nGs = size(mpc.Gsou,1);
nGb = size(mpc.Gbus,1);
nGl = size(mpc.Gline,1);
nPTG = size(mpc.ptg,1);

%
% PGs = PGs * 1e6 / 24 / 3600; % m3/s
% Pg = Pg * baseMVA * 1e6; % w
% Qptg = Qptg * 1e6/24/3600;


%%
for r = 1:nGasType
    % for each type of gas
    PGsbus = mpc.Gsou(:,1) ; 
    Cgs_PGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
    Qdbus = iGasLoad; 
    Cgs_Qd = sparse(Qdbus, (1:nGasLoad)', 1, nGb, nGasLoad); % connection matrix
    f_r = Cgs_PGs*(PGs .* mpc.gasCompositionForGasSource(:,r)) - Cgs_Qd * Qd(:,r); % supply-demand, Mm3/day
    
    % gas flow
    for  m = 1:nGl
        fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);

        f_r(fb) = f_r(fb) - gasFlow(m,r);
        f_r(tb) = f_r(tb) + gasFlow(m,r);
    end
    % ptg
    if (r == 1) || (r == 5) % is methane or hydrogen
        for i = 1:nPTG
            GB = mpc.ptg(i,1);
            if r == 1 % methane
                ptgGasType = 1;
            elseif r == 5 % hydrogen
                ptgGasType = 2;
            end
            f_r(GB) = f_r(GB) + Qptg(i,ptgGasType);
        end
    end
    % gfu
    for i = 1:nGPP
        GB = mpc.GEcon(find(mpc.GEcon(:,2) == mpc.gen(mpc.gfuIndex(i),1)),1);
        f_r(GB) = f_r(GB)-Qgpp(i,r);
    end
    f(:,r) = f_r;
end

end