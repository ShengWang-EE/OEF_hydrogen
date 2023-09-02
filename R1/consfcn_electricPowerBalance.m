function g = consfcn_electricPowerBalance(Va,Pg,Pptg,mpc)
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
nb = length(Va);             %% number of buses
ng = length(Pg);            %% number of dispatchable injections
nPTG = size(mpc.ptg,1);
%% unpack data
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
%%
% gen(:, PG) = Pg * baseMVA;  %% active generation in MW
  [B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
  neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   %% Pbus w.r.t. Pg
  Amis = [B neg_Cg];
  %-----modified-------
  bmis = -(bus(:, PD) + bus(:, GS)) / baseMVA - Pbusinj;
  g = Amis*[Va;Pg]-bmis;%MW/baseMVA
  % ptg
  for i = 1:nPTG
      EBindex = mpc.ptg(i,2);
      g(EBindex) = g(EBindex) + Pptg(i); %待定，应该是负荷为正？
  end

end