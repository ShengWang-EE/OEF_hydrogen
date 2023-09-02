function h = consfcn_electricBranchFlow(Va, mpc, il)
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
    upf = mpc.branch(il, RATE_A) / mpc.baseMVA - Pfinj(il);
    upt = mpc.branch(il, RATE_A) / mpc.baseMVA + Pfinj(il);
%% unpack data
h1 = Bf(il,:)*Va - upf;
h2 = -upt - Bf(il,:)*Va;
h = [h1;h2];
end