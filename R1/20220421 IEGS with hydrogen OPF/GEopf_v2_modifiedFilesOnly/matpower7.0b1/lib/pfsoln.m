function [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, Ybus, Yf, Yt, V, ref, pv, pq, mpopt)
%PFSOLN  Updates bus, gen, branch data structures to match power flow soln.
%   [BUS, GEN, BRANCH] = PFSOLN(BASEMVA, BUS0, GEN0, BRANCH0, ...
%                                   YBUS, YF, YT, V, REF, PV, PQ, MPOPT)

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% default options
if nargin < 12
    mpopt = mpoption();
end

%% initialize return values
bus     = bus0;
gen     = gen0;
branch  = branch0;

%%----- update bus voltages -----
bus(:, VM) = abs(V);
bus(:, VA) = angle(V) * 180 / pi;

%%----- update Qg for gens at PV/slack buses and Pg for slack bus(es) -----
%% generator info
on = find(gen(:, GEN_STATUS) > 0 & ...  %% which generators are on?
        bus(gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and not at PQ buses
off = find(gen(:, GEN_STATUS) <= 0);    %% which generators are off?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

%% compute total injected bus powers
Sbus = V(gbus) .* conj(Ybus(gbus, :) * V);

%% update Qg for generators at PV/slack buses
gen(off, QG) = zeros(length(off), 1);   %% zero out off-line Qg
%% don't touch the ones at PQ buses
[Pd_gbus, Qd_gbus] = total_load(bus(gbus, :), [], 'bus', [], mpopt);
gen(on, QG) = imag(Sbus) * baseMVA + Qd_gbus;   %% inj Q + local Qd
%% ... at this point any buses with more than one generator will have
%% the total Q dispatch for the bus assigned to each generator. This
%% must be split between them. We do it first equally, then in proportion
%% to the reactive range of the generator.

if length(on) > 1
    %% build connection matrix, element i, j is 1 if gen on(i) at bus j is ON
    nb = size(bus, 1);
    ngon = size(on, 1);
    Cg = sparse((1:ngon)', gbus, ones(ngon, 1), ngon, nb);

    %% divide Qg by number of generators at the bus to distribute equally
    ngg = Cg * sum(Cg)';    %% ngon x 1, number of gens at this gen's bus
    gen(on, QG) = gen(on, QG) ./ ngg;

    %% prep 
    Qmin = gen(on, QMIN);
    Qmax = gen(on, QMAX);
    M = abs(gen(on, QG));
    M(~isinf(Qmax)) = M(~isinf(Qmax)) + abs(Qmax(~isinf(Qmax)));
    M(~isinf(Qmin)) = M(~isinf(Qmin)) + abs(Qmin(~isinf(Qmin)));
    M = Cg * Cg' * M;
    Qmin(Qmin ==  Inf) =  M(Qmin ==  Inf);
    Qmin(Qmin == -Inf) = -M(Qmin == -Inf);
    Qmax(Qmax ==  Inf) =  M(Qmax ==  Inf);
    Qmax(Qmax == -Inf) = -M(Qmax == -Inf);

    %% divide proportionally
    Cmin = sparse((1:ngon)', gbus, Qmin, ngon, nb);
    Cmax = sparse((1:ngon)', gbus, Qmax, ngon, nb);
    Qg_tot = Cg' * gen(on, QG);     %% nb x 1 vector of total Qg at each bus
    Qg_min = sum(Cmin)';            %% nb x 1 vector of min total Qg at each bus
    Qg_max = sum(Cmax)';            %% nb x 1 vector of max total Qg at each bus
    gen(on, QG) = Qmin + ...
        (Cg * ((Qg_tot - Qg_min)./(Qg_max - Qg_min + eps))) .* ...
            (Qmax - Qmin);          %%                ^ avoid div by 0

    %% fix gens at buses with Qg range = 0 (use equal violation for all)
    ig = find(abs(Cg * (Qg_min - Qg_max)) < 10*eps);  %% gens at buses with Qg range = 0
    if ~isempty(ig)
        ib = find(sum(Cg(ig,:), 1)');   %% buses with Qg range = 0
        %% total mismatch @ bus div by number of gens
        mis = sparse(ib, 1, (Qg_tot(ib) - Qg_min(ib)) ./ sum(Cg(:, ib)', 2), nb, 1);
        gen(on(ig), QG) = Qmin(ig) + Cg(ig, :) * mis;
    end
end                                             %% (terms are mult by 0 anyway)

%% update Pg for slack gen(s)
for k = 1:length(ref)
    refgen = find(gbus == ref(k));              %% which is(are) the reference gen(s)?
    Pd_refk = total_load(bus(ref(k), :), [], 'bus', [], mpopt);
    gen(on(refgen(1)), PG) = real(Sbus(refgen(1))) * baseMVA + Pd_refk; %% inj P + local Pd
    if length(refgen) > 1       %% more than one generator at this ref bus
        %% subtract off what is generated by other gens at this bus
        gen(on(refgen(1)), PG) = gen(on(refgen(1)), PG) ...
                                - sum(gen(on(refgen(2:length(refgen))), PG));
    end
end

%%----- update/compute branch power flows -----
out = find(branch(:, BR_STATUS) == 0);      %% out-of-service branches
br = find(branch(:, BR_STATUS));            %% in-service branches
Sf = V(branch(br, F_BUS)) .* conj(Yf(br, :) * V) * baseMVA; %% complex power at "from" bus
St = V(branch(br, T_BUS)) .* conj(Yt(br, :) * V) * baseMVA; %% complex power injected at "to" bus
branch(br, [PF, QF, PT, QT]) = [real(Sf) imag(Sf) real(St) imag(St)];
branch(out, [PF, QF, PT, QT]) = zeros(length(out), 4);
