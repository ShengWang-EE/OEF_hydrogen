function f = consfcn_gasBalance(PGs,Pg,composition_hy,composition_gas,Qptg, gasFlow, signGf,mpc, GCV_hy, GCV_gas,etaGFU,baseMVA)
%% parameter
nGs = size(mpc.Gsou,1);
nGb = size(mpc.Gbus,1);
nGl = size(mpc.Gline,1);
nPTG = size(mpc.ptg,1);
nGFU = size(mpc.gfuIndex,1);

%
PGs = PGs * 1e6 / 24 / 3600; % m3/s
Pg = Pg * baseMVA * 1e6; % w
Qptg = Qptg * 1e6/24/3600;
energyProduction = PGs * GCV_gas; %w
energyDemand = mpc.Gbus(:,3) * 1e6/24/3600 * GCV_gas;

%%
PGsbus = mpc.Gsou(:,1); 
Cgs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
f = Cgs*energyProduction - energyDemand; % supply-demand

% gas flow
for  m = 1:nGl
    if signGf(m) == 1 % positive
        fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);
    else
        fb = mpc.Gline(m,2); tb = mpc.Gline(m,1);
    end
    GCV = GCV_hy * composition_hy(fb) + GCV_gas*composition_gas(fb);
    f(fb) = f(fb) - gasFlow(m) * GCV;
    f(tb) = f(tb) + gasFlow(m) * GCV;
end
% ptg
for i = 1:nPTG
    GB = mpc.ptg(i,1);
    f(GB) = f(GB) + Qptg(i) * GCV_hy;
end
% gfu
for i = 1:nGFU
    GB = mpc.GEcon(find(mpc.GEcon(:,2) == mpc.gen(mpc.gfuIndex(i),1)),1);
    f(GB) = f(GB)-Pg(mpc.gfuIndex(i)) / etaGFU;
end



end