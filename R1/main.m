clear
clc
yalmip('clear')
%% parameters
[mpc,gtd] = case24GEv7(); 
mpc.gen(:,10) = 0; % min output zero

%
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
baseMVA = 100;
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
%
nb   = size(mpc.bus, 1);                                                    %% number of buses
nGb  = size(mpc.Gbus,1);                                                    % number of gas bus
nGl = size(mpc.Gline,1);
nGen = sum(mpc.gen(:,22)==1)+sum(mpc.gen(:,22)==0);                         % all gen(TFU and GFU), excluded dispatchable loads
nGPP = size(mpc.gfuIndex,1);
nGs = size(mpc.Gsou,1);
nGasLoad = size(find(mpc.Gbus(:,3)~=0),1);
nPTG = size(mpc.ptg,1);
iGasLoad = find(mpc.Gbus(:,3)~=0);

% CH4, C2H6, C3H8, C4H10, H2, N2, CO2
nGasType = 7; iCombustibleGas = 1:5; iNonCombustibleGas = 6:7;
[GCVall, Mall, M_air, fsAll, aAll, R_air, T_stp, Prs_stp, Z_ng, ...
    T_gas, eta_electrolysis,eta_methanation,etaGFU] = initializeParameters();
%
xi = 0.1;                                                                   % threshold for the security index
alpha_PHI = 1;alpha_x = 1; alpha_Qd = 1;
lambda = 10;                                                                % coefficient for gas flow
multiplier = 2; alpha_PHI_max = 1e4; alpha_x_max = 1e4; alpha_Qd_max = 1e4;
iterationMax = 40;
gap = 1e-2;                                                                 % optimization convergence criterion

ob.Pg = [];

% step2: get initial results
GEresult0 = GErunopf(mpc);

%% initial value
Prs_square0 = GEresult0.Gbus(:,7).^2;
PGs0 = GEresult0.Gsou(:,5);
Qptg0 = 0;
Pg0 = GEresult0.gen(:,2)/baseMVA;
Va0 = GEresult0.bus(:,9)/180*pi;
LCg0 = 0;

% parameters for original natural gas
gasComposition_ng = mean(mpc.gasCompositionForGasSource);
M_ng = gasComposition_ng * Mall';
GCV_ng = gasComposition_ng * GCVall';
S_ng = M_ng/M_air;
WI_ng = GCV_ng / sqrt(S_ng);
FS_ng = gasComposition_ng * fsAll';
CP_ng = gasComposition_ng * aAll' / sqrt(S_ng);

% variables that need precalculation and iteration
S_pipeline0 = M_ng/M_air; %
Z0 = Z_ng;
gasComposition0 = repmat(gasComposition_ng,[nGb,1]);
gasFlow_sum0 = GEresult0.Gline(:,6);
QptgInbus0 = zeros(nGb,1);
signGf = sign(GEresult0.Gline(:,6));
W0 = nodalGasInjection(PGs0,QptgInbus0,signGf, gasFlow_sum0,mpc,nGs,nGb,nGl);
gasLoad_ng = mpc.Gbus(iGasLoad,3);
%% case settings
% general, apply to all cases
mpc.gencost(1:4,2:7) = 0;
mpc.branch([28,31,32,33],6) = mpc.branch([28,31,32,33],6)/5;                % line capacity surround RNG / 5 围绕着RNG的线路容量/5
subsidyFlag = 1;
% alpha_PHI = 1e-2; alpha_x = 1e-2; alpha_Qd = 1e-2;
% alpha_PHI_max = 1e3; alpha_x_max = 1e4; alpha_Qd_max = 1e4;
% --- base case: Section VI.A and VI.B， connected gas bus for ptg is 1,4,10, the original case 认为ptg连接天然气节点的编号是1，4，10，就是原case不用修改
% --- comment 2-1: ---
% xi = 0.3;
% mpc.ptg(:,6) = 0; % case A: no hydrogen
% mpc.ptg(:,6) = 0.5; % case B
% mpc.ptg(:,6) = 2; % case C
% mpc.ptg(:,6) = 5; % case D

% --- comment 2-3: add adjustive gas ---
% adjustiveGas.Gsou = [
%     1    0   0   2
%     1    0   0   2
%     4    0   0   2
%     4    0   0   2
%     10   0   0   2
%     10   0   0   2
% ];
% adjustiveGas.gasCompositionForGasSource = [
%     91.1	4.3	3	1.4	0	0.2	    0
%     0	0	0	0	0	    100	    0
%     91.1	4.3	3	1.4	0	0.2	    0
%     0	0	0	0	0	    100	    0
%     91.1	4.3	3	1.4	0	0.2	    0
%     0	0	0	0	0	    100	    0
% ]/100;
% mpc.Gsou = [mpc.Gsou; adjustiveGas.Gsou];
% mpc.gasCompositionForGasSource = [mpc.gasCompositionForGasSource; adjustiveGas.gasCompositionForGasSource];
% nGs = size(mpc.Gsou,1);
% alpha_PHI_max = 1e2; alpha_x_max = 1e2; alpha_Qd_max = 1e2;
% case A:
% xi = 0.05;  mpc.ptg(:,6) = 0.5;   priceFactor = 1.5;
% case B
% xi = 0.05;  mpc.ptg(:,6) = 2;   priceFactor = 1.5;
% case C
% xi = 0.10;  mpc.ptg(:,6) = 2;   priceFactor = 1.5;
% case D
% xi = 0.05;  mpc.ptg(:,6) = 2;   priceFactor = 3;
% adjustiveGas.Gcost = [
%     3541.66666666667*priceFactor; 3541.66666666667*priceFactor; 3541.66666666667*priceFactor;
%     3541.66666666667*priceFactor; 3541.66666666667*priceFactor; 3541.66666666667*priceFactor;
%     ];
% mpc.Gcost = [mpc.Gcost; adjustiveGas.Gcost;];

% --- comment 3-3: PTG subsidy ---
% mpc.ptg = [
%     1	10	2	200	0	2
%     8	17	2	200	0	2
%     5	18	2	200	0	2
%     ];
% case A
% subsidyFlag = 1;
% case B
% subsidyFlag = 0; 
% case C 
% subsidyFlag = 0; 
% mpc.gencost(mpc.windfarmIndex,:) = 	[2	1500	0	3	0	4.4231	395.3749];
% mpc.branch([28,31,32,33],6) = mpc.branch([28,31,32,33],6) *5;

% --- comment 4-14: impact of uncertainties ---
% mpc.ptg(:,6) = 3;
% case A
% mpc.gen(mpc.windfarmIndex,PMAX) = mpc.gen(mpc.windfarmIndex,PMAX) * 0.5; % with security cons
% mpc.gen(mpc.windfarmIndex,PMAX) = mpc.gen(mpc.windfarmIndex,PMAX) * 0.5; xi = 1;% without security cons
% case B
% mpc.gen(mpc.windfarmIndex,PMAX) = mpc.gen(mpc.windfarmIndex,PMAX) * 1; % with security cons
% mpc.gen(mpc.windfarmIndex,PMAX) = mpc.gen(mpc.windfarmIndex,PMAX) * 1; xi = 1;% without security cons
% case C
% mpc.gen(mpc.windfarmIndex,PMAX) = mpc.gen(mpc.windfarmIndex,PMAX) * 2; % with security cons
% mpc.gen(mpc.windfarmIndex,PMAX) = mpc.gen(mpc.windfarmIndex,PMAX) * 2; xi = 1;% without security cons
% case D
% mpc.Gbus(:,3) = mpc.Gbus(:,3) * 0.8;
% case E
% mpc.Gbus(:,3) = mpc.Gbus(:,3) * 1;
% case F
% mpc.Gbus(:,3) = mpc.Gbus(:,3) * 1.05;

% --- comment 4-15 ---
mpc.ptg(:,6) = 3;
% case A: base case
% case B:
% mpc.branch(30,RATE_A:RATE_C) = 0;
% case C 

% case D: same as case A
% case E:
% mpc.ptg = [
%     16	10	2	200	0	3
%     20	17	2	200	0	3
%     19	18	2	200	0	3
%     ];
% case F
% mpc.gen(mpc.windfarmIndex,1) = 1;
% case G:
mpc.GEcon=[
    15 15
    16 1
    19 7
    20 2
];


% mpc.gen(23,9) = 200;
% mpc.gen(1:4,9) = 50;
% mpc.ptg = [
%     1	10	2	200	0	0.5
%     8	17	2	200	0	0.5
%     5	18	2	200	0	0.5
%     ];
%%
solverTime = zeros(iterationMax,1);
for v = 1:iterationMax
%% state variables
Prs_square = sdpvar(nGb,1); % bar^2
PGs = sdpvar(nGs,1); % Mm3/day
Qd = sdpvar(nGasLoad,nGasType);% Mm3/day
Qptg = sdpvar(nPTG,2); % [ methane; hydrogen ] % Mm3/day
Pptg = sdpvar(nPTG,1); % electricity consumption, 1/100 MW
Pg = sdpvar(size(mpc.gen,1),1); % include TPP, GPP and renewable generators, 1/100 MW
Qgpp = sdpvar(nGPP, nGasType); % Mm3/day
Va = sdpvar(nb,1);
gamma = sdpvar(nGl,1); % direction of the gas flow
gasComposition = sdpvar(nGb,nGasType); 
gasFlow = sdpvar(nGl,nGasType);% Mm3/day
PHI = sdpvar(nGl,1); %auxiliary variable for gas flow
sigma_PHI = sdpvar(nGl,1); % error limit for gas flow
varepsilon_x = sdpvar(nGb, nGasType); % error for gas composition
sigma_x1000 = sdpvar(nGb,nGasType); % error limit for gas composition
sigma_Qd = sdpvar(nGasLoad, nGasType);
% ------------test------------------
LCg = zeros(nGasLoad,1);
% varepsilon_x = zeros(nGb, nGasType); % error for gas composition
gamma = signGf;
%% bounds
Prsmin = mpc.Gbus(:,5); Prsmax = mpc.Gbus(:,6); % bar
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4); % Mm3/day
Qptgmin = mpc.ptg(:,5); Qptgmax = mpc.ptg(:,6);
QptgMax_hydrogen = Qptgmax;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
Pgmax(34:end,:) = 1e-4; % no electricity load curtailment
LCgmin = zeros(nGasLoad,1);
LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*0;  
refs = find(mpc.bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = 1;   %% voltage angle reference constraints
Val(refs) = 1;
gasFlowMax = mpc.Gline(:,5);
%% pre-calculation
% calculate the gas flow of each pipeline
C = mpc.Gline(:,3) .* sqrt(S_ng)*sqrt(Z_ng) ./ sqrt(S_pipeline0) ./ sqrt(Z0);


%% constraints
% gas demand cons
energyDemand = mpc.Gbus(iGasLoad,3) * GCV_ng; % energy need of these gas bus
QdGasCompositionCons = [];
for ii = 1:nGasLoad
    QdGasCompositionCons = [
        QdGasCompositionCons;
        (gasComposition0(ii,:)) .* sum(Qd(ii,:),2) - sigma_Qd(ii,:) <= Qd(ii,:) ...
            <= (gasComposition0(ii,:)) .* sum(Qd(ii,:),2) + sigma_Qd(ii,:); % reformulated
%         gasComposition(ii,:) .* sum(Qd(ii,:),2) == Qd(ii,:); % original cons
        ];
end
gasDemandCons = [
    QdGasCompositionCons;
    Qd * GCVall'/1e9 == energyDemand/1e9;    
    Qd >= 0;
    0 <= sigma_Qd <= repmat(gasLoad_ng,[1,nGasType]);
    
%     PGsmin <= PGs <= PGsmax;
    ]:'gasDemandCons';
% nodal gas flow balance cons
nodalGasFlowBalanceCons = [
    consfcn_nodalGasFlowBalance(PGs,Qd,Qgpp,Qptg, gasFlow,mpc,nGasType,nGPP,nGasLoad,iGasLoad) == 0;
    ]:'nodalGasFlowBalanceCons';
% ptg
PTGcons = [
    ( Qptg(:,1) * 1e6/24/3600 *GCVall(1) / eta_methanation + Qptg(:,2) * 1e6/24/3600 * GCVall(5) ...
        ) /1e6 == Pptg*baseMVA / eta_electrolysis; % w
    0 <= Pptg*baseMVA / eta_electrolysis <= QptgMax_hydrogen /24/3600 * GCVall(5); 
    0 <= Qptg;
    ];
% gpp
Pgpp = Pg(mpc.gfuIndex) * baseMVA;% MW
GPPcons = [
    Pgpp == Qgpp/24/3600 * GCVall';
    Qgpp >= 0;
    ]:'GPPcons';
% electricity flow
electricityCons = [
    consfcn_electricPowerBalance(Va,Pg,Pptg,mpc) == 0;
    consfcn_electricBranchFlow(Va, mpc, il) <= 0;
    Pgmin <= Pg <= Pgmax;
    Va(refs) == 0; 
    ]:'electricityCons';
% wobbe index
GCV_nodal = gasComposition * GCVall';
S_nodal = gasComposition * Mall' / M_air;
sqrtS = 0.5 * (S_nodal/sqrt(S_ng) + sqrt(S_ng));
WI_nodal = GCV_nodal ./ sqrtS;                                              % if can't be auto converted, then manually converted如果不能自动转化，那就手动化一下
WobbeIndexCons = [
    (1-xi) * WI_ng * sqrtS/1e6 <= GCV_nodal/1e6 <= (1+xi) * WI_ng * sqrtS/1e6
    ]:'WobbeIndexCons';
% Weaver flame speed factor
FSnodal = gasComposition * fsAll';
FScons = [
    (1-xi) * FS_ng <= FSnodal <= (1+xi) * FS_ng;
    ]:'FScons';
% combustion potential
CPnodal = gasComposition * aAll' ./ sqrtS; 
CPcons = [
    (1-xi) * CP_ng * sqrtS <= gasComposition * aAll' <= (1+xi) * CP_ng * sqrtS;
    ]:'CPcons';
% other security cons
otherSecurityCons = [
    0 <= gasComposition <= 1;
    gasComposition(:,5) <= xi;
    (1-xi)*GCV_ng/1e6 <= GCV_nodal/1e6 <= (1+xi)*GCV_ng/1e6;
    ]:'otherSecurityCons';
% SOC reformulation for gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
gasFlow_sum = sum(gasFlow,2);
gasFlowSOCcons = [
    (gamma-1) .* gasFlowMax / 2 <= gasFlow_sum <= (gamma+1) .* gasFlowMax / 2;
    repmat((gamma-1) .* gasFlowMax / 2,[1,nGasType]) <= gasFlow <= repmat((gamma+1) .* gasFlowMax / 2, [1,nGasType]);
    Prsmin.^2 <= Prs_square <= Prsmax.^2;
    PHI >= gasFlow_sum.^2 ./ C.^2 ;
    PHI <= gasFlow_sum0 + 2*gasFlow_sum0 .* (gasFlow_sum - gasFlow_sum0) ./ C.^2 + sigma_PHI;
%     PHI == gasFlow_sum.* gasFlow_sum0./ C.^2; % original cons
    PHI >= Prs_square(TB) - Prs_square(FB) + (gamma + 1) .* (Prsmin(FB).^2 - Prsmax(TB).^2);
    PHI >= Prs_square(FB) - Prs_square(TB) + (gamma - 1) .* (Prsmax(FB).^2 - Prsmin(TB).^2);
    PHI <= Prs_square(TB) - Prs_square(FB) + (gamma + 1) .* (Prsmax(FB).^2 - Prsmin(TB).^2);
    PHI <= Prs_square(FB) - Prs_square(TB) + (gamma - 1) .* (Prsmin(FB).^2 - Prsmax(TB).^2);

    1e4 >= sigma_PHI >= 0;
%     PHI <= 1e4;
    ]:'gasFlowSOCcons';
% gas composition Taylor
PTGbus = mpc.ptg(:,1) ; 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbusMethane = CgsPTG * Qptg(:,1) ; QptgInbusHydrogen = CgsPTG * Qptg(:,2);
QptgInbusForAllGasComposition = [QptgInbusMethane, zeros(nGb,3),QptgInbusHydrogen,zeros(nGb,2)];
gasCompositionCons1 = [];
nodalGasInjectionSum = 0;
for r = 1:nGasType
    nodalGasInjectionSum = nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
            gamma, gasFlow(:,r),mpc,nGs,nGb,nGl);
end
for r = 1:nGasType
    gasCompositionCons1 = [
        gasCompositionCons1;
        gasComposition(:,r) == ...
            nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
            gamma, gasFlow(:,r),mpc,nGs,nGb,nGl) ./ W0 + varepsilon_x(:,r);
%         gasComposition(:,r) == nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),QptgInbusForAllGasComposition(:,r),...
%             gamma, gasFlow(:,r),mpc,nGs,nGb,nGl) ./ nodalGasInjectionSum; % original cons
        ];
end
gasCompositionCons = [
    gasCompositionCons1;
    -sigma_x1000/1000 <= varepsilon_x <= sigma_x1000/1000;
    1000 >= sigma_x1000 >= 0;
    sum(gasComposition,2) == 1;
    gasComposition(12,:) == gasComposition(17,:);
    ]:'gasCompositionCons';
% summarize all the cons
constraints = [
    gasDemandCons;
    nodalGasFlowBalanceCons;
    PTGcons;
    GPPcons;
    electricityCons;
    WobbeIndexCons;
    FScons;
    CPcons;
    otherSecurityCons;
    gasFlowSOCcons;
    gasCompositionCons;
    PGsmin <= PGs <= PGsmax;
    ];
%% solve the problem
objfcn = obj_operatingCost(Pg,PGs,LCg,Qptg,PHI, sigma_PHI, sigma_x1000, alpha_PHI,alpha_x, mpc,lambda, subsidyFlag) ...
    + alpha_Qd * 100 * sum(sum(sigma_Qd));
options = sdpsettings('verbose',2,'solver','gurobi', 'debug',1,'usex0',0);
% options.ipopt.tol = 1e-4;
output{v} = optimize(constraints, objfcn, options);

solverTime(v) = output{v}.solvertime;
%% results
Prs_square = value(Prs_square);
Prs = sqrt(value(Prs_square));
PGs = value(PGs); % Mm3/day
Qd = value(Qd); % Mm3/day
Qptg = value(Qptg);
Pptg = value(Pptg) * baseMVA; % MW
Pg = value(Pg) * baseMVA;  % MW
Pgpp = value(Pgpp) * baseMVA;  % MW
Qgpp = value(Qgpp);  % Mm3/day
Va = value(Va);
LCg = value(LCg);
gamma = value(gamma);
gasComposition = value(gasComposition); % hydrogen, gas
gasFlow = value(gasFlow);
gasFlow_sum = sum(value(gasFlow_sum),2);
S_nodal = value(S_nodal);
GCV_nodal = value(GCV_nodal)/1e6; % MJ/m3
WInodal = value(WI_nodal)/1e6; % MJ/m3
FSnodal = value(FSnodal);
CPnodal = value(CPnodal);
PHI = value(PHI);
sigma_PHI = value(sigma_PHI);
varepsilon_x = value(varepsilon_x);
sigma_x1000 = value(sigma_x1000);
sigma_x = value(sigma_x1000)/1000;
sigma_Qd = value(sigma_Qd);

[objfcn,totalCost,genAndLCeCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy,penalty_PHI,penalty_sinma_PHI, penalty_sigma_x1000] = ...
    obj_operatingCost(Pg,PGs,LCg,Qptg,PHI, sigma_PHI, sigma_x1000, alpha_PHI,alpha_x, mpc,lambda, subsidyFlag);
%
[sol{v}.totalCost,sol{v}.genAndLCeCost,sol{v}.gasPurchasingCost,sol{v}.gasCurtailmentCost,sol{v}.PTGsubsidy,sol{v}.penalty_PHI,sol{v}.penalty_sinma_PHI, sol{v}.penalty_sigma_x, sol{v}.objfcn, ...
    sol{v}.Prs_square,sol{v}.Prs,sol{v}.PGs,sol{v}.Qd,sol{v}.Qptg,sol{v}.Pptg,sol{v}.Pg,sol{v}.Pgpp,sol{v}.Qgpp,...
    sol{v}.Va,sol{v}.LCg,sol{v}.gamma,sol{v}.gasComposition,sol{v}.gasFlow,sol{v}.gasFlow_sum,sol{v}.S_nodal,...
    sol{v}.GCV_nodal,sol{v}.WI,sol{v}.PHI,sol{v}.sigma_PHI,sol{v}.varepsilon_x,sol{v}.sigma_x,sol{v}.sigma_Qd] = ...
    deal(totalCost,genAndLCeCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy,penalty_PHI,penalty_sinma_PHI, penalty_sigma_x1000, objfcn, ...
    Prs_square,Prs,PGs,Qd,Qptg,Pptg,Pg,Pgpp,Qgpp,Va,LCg,gamma,gasComposition,gasFlow,gasFlow_sum,S_nodal,...
    GCV_nodal,WInodal,PHI,sigma_PHI,varepsilon_x,sigma_x,sigma_Qd);

%% covergence criterion
if v > 1
    criterion.PHI(v-1) = abs( ( sum(sum(sigma_PHI)) - sum(sum(sol{v-1}.sigma_PHI)) ) / sum(sum(sol{v-1}.sigma_PHI)) );
    criterion.x(v-1) = abs( ( sum(sum(sigma_x)) - sum(sum(sol{v-1}.sigma_x)) ) / sum(sum(sol{v-1}.sigma_x)) );
    %cost
    criterion.totalCost(v-1) = abs( (totalCost - sol{v-1}.totalCost)/sol{v-1}.totalCost );
    criterion.genAndLCeCost(v-1) = abs( (genAndLCeCost - sol{v-1}.genAndLCeCost)/sol{v-1}.genAndLCeCost );
    criterion.gasPurchasingCost(v-1) = abs( (gasPurchasingCost - sol{v-1}.gasPurchasingCost)/sol{v-1}.gasPurchasingCost );
    criterion.PTGsubsidy(v-1) = abs( (PTGsubsidy - sol{v-1}.PTGsubsidy)/sol{v-1}.PTGsubsidy );

    criterion.S_nodal(v-1) = max(abs( (S_nodal-sol{v-1}.S_nodal)./sol{v-1}.S_nodal ));
    criterion.Z(v-1) = 0;
    criterion.Qd(v-1) = abs( ( sum(sum(sigma_Qd)) - sum(sum(sol{v-1}.sigma_Qd)) ) / sum(sum(sol{v-1}.sigma_Qd)) );

    if (criterion.PHI(v-1)<=gap) && (criterion.x(v-1)<=gap) && (criterion.totalCost(v-1)<=gap) ...
            && (criterion.S_nodal(v-1)<=gap) && (criterion.Z(v-1)<=gap) && (criterion.Qd(v-1)<=gap) 
        break
    else
        alpha_PHI = min([multiplier * alpha_PHI,alpha_PHI_max]);
        alpha_x = min([multiplier * alpha_x,alpha_x_max]);
        alpha_Qd = min([multiplier * alpha_Qd,alpha_Qd_max]);
    end
end
% ob
ob.Pg = [ob.Pg, Pg];
%% update S and Z

S_pipeline0 = ( (1+gamma).*S_nodal(FB) + (1-gamma).*S_nodal(TB) ) / 2;
Z0 = Z_ng; % consider Z is the same
QptgInbus = CgsPTG * sum(Qptg,2);
W0 = nodalGasInjection(PGs,QptgInbus,gamma,gasFlow_sum,mpc,nGs,nGb,nGl);
for r = 1:nGasType
    gasComposition0(:,r) = nodalGasInjection(PGs.*mpc.gasCompositionForGasSource(:,r),...
        QptgInbusForAllGasComposition(:,r), gamma, gasFlow(:,r),mpc,nGs,nGb,nGl) ./ W0 + varepsilon_x(:,r);
end
if output{v}.problem ~= 0
%     error('optimization failed');
%     break
end
end
securityIndices = [WInodal,FSnodal,CPnodal,GCV_nodal];
save






