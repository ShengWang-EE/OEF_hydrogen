clear
clc
yalmip('clear')
% 经验：先确定潮流方向的确很重要，能快很多
%%
[GCVall, Mall, M_air, fsAll, aAll, R_air, T_stp, Prs_stp, Z_ng, ...
    T_gas, eta_electrolysis,eta_methanation,etaGFU] = initializeParameters();
mpc = caseUKgasSystem();
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl = size(mpc.Gline,1);
nGs = size(mpc.Gsou,1);
nGasLoad = size(find(mpc.Gbus(:,2)~=0),1);
iGasLoad = find(mpc.Gbus(:,2)~=0);

% CH4, C2H6, C3H8, C4H10, H2, N2, CO2
nGasType = 7; 
%%
% var
PGs = sdpvar(nGs,1);
gasflow = sdpvar(nGl,nGasType);
% gamma = binvar(nGl,1);
gamma = mpc.gamma;
Qd = sdpvar(nGasLoad,nGasType);% Mm3/day
gasComposition = sdpvar(nGb,nGasType); 
LCg = sdpvar(nGasLoad,1);
slack_Qd = sdpvar(nGasLoad,nGasType); 
slack_nodal = sdpvar(nGb,nGasType); 
auxiliaryVar = sdpvar(nGl,nGasType);
%

flowDirection = (gamma-0.5) * 2; % {0,1}->{-1,+1}
PGsmin = 0; PGsmax = mpc.Gsou(:,2)/2;
gasflowmax = 400;
energyDemand = mpc.Gbus(iGasLoad,2);
gasComposition0 = sum(repmat(mpc.Gsou(:,2),[1,nGasType]) .* mpc.gasComposition,1)./sum(sum(repmat(mpc.Gsou(:,2),[1,nGasType]) .* mpc.gasComposition,1));
penalty_Qd = 100;

%% cons
boxCons = [
    PGsmin <= PGs <= PGsmax;
    gasflowmax >= gasflow >= 0; % 全是正，gamma>0则是顺应pipe方向，反之则反过来
    Qd >= 0;
    0 <= gasComposition <= 1;
    sum(gasComposition,2) == 1;
    energyDemand >= LCg >= 0;
    gasComposition(:,5) <= 0.1;
    0 <= gamma <= 1;
    ];
% gas demand cons
QdCons = [
    Qd == gasComposition(iGasLoad,:) .* repmat(sum(Qd,2),[1,nGasType]);
%     -slack_Qd <= repmat(gasComposition0,[nGasLoad,1]) .* repmat(sum(Qd,2),[1,nGasType]) - Qd <= slack_Qd;
    slack_Qd >= 0;   
    Qd * GCVall'/1e6 + LCg >= energyDemand; % MJ
    ];

% nodal gas balance
for r = 1:nGasType
    PGsbus = mpc.Gsou(:,1) ; 
    Cgs_PGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
    Qdbus = iGasLoad; 
    Cgs_Qd = sparse(Qdbus, (1:nGasLoad)', 1, nGb, nGasLoad); % connection matrix
    f_r(:,r) = Cgs_PGs*(PGs .* mpc.gasComposition(:,r)) - Cgs_Qd * Qd(:,r); % supply-demand, Mm3/day    
    % gas flow
    for  m = 1:nGl
        fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);
        f_r(fb,r) = f_r(fb,r) - flowDirection(m) * gasflow(m,r);
        f_r(tb,r) = f_r(tb,r) + flowDirection(m) * gasflow(m,r);
    end
end
nodalGasBalanceCons = [f_r == 0;];
% gas composition cons
for i = 1:nGb
    asFromBus_PipeIndex = find(mpc.Gline(:,1) == i); asToBus_PipeIndex = find(mpc.Gline(:,2) == i);
    assumingfrombus = mpc.Gline(asToBus_PipeIndex,1); assumingtobus = mpc.Gline(asFromBus_PipeIndex,2);
%     nodalgasInjection(i,:) = sum(gasflow(asFromBus_PipeIndex,:) .* repmat(gamma(asFromBus_PipeIndex),[1,nGasType]),1) ...
%         + sum(gasflow(asToBus_PipeIndex,:) .* (1-repmat(gamma(asToBus_PipeIndex),[1,nGasType])),1);
    nodalgasInjection(i,:) = sum(auxiliaryVar(asFromBus_PipeIndex,:),1) ...
        + sum(gasflow(asToBus_PipeIndex,:) - auxiliaryVar(asToBus_PipeIndex,:),1);
end
nodalGasCompositionCons = [
    auxiliaryVar <= gasflow;
    auxiliaryVar >= gasflow - (1-repmat(gamma,[1,nGasType])).* gasflowmax;
    gasflowmax >= auxiliaryVar >= 0;
    nodalgasInjection == gasComposition .* repmat(sum(nodalgasInjection,2),[1,nGasType]);
%     -slack_nodal <= nodalgasInjection - repmat(gasComposition0,[nGb,1]) .* repmat(sum(nodalgasInjection,2),[1,nGasType]) <= slack_nodal;
    slack_nodal >= 0;
];
% 所有下游管道的composition等于直接上游节点的gascomposition
% for m = 1:nGl
%     gasCompositionPipe = gasflow(m,:)./ repmat(sum(gasflow(m,:),2),[1,nGasType]);
%     gasCompositionPipe == mpc.Gline(m,1) * gamma(m) + mpc.Gline(m,2) * (1-gamma(m));
% end

pipeGasCompositionCons = [
    gasflow == ( gasComposition(mpc.Gline(:,1),:) .* repmat(gamma,[1,nGasType]) + gasComposition(mpc.Gline(:,2),:) .* (1-repmat(gamma,[1,nGasType])) ) .* repmat(sum(gasflow,2),[1,nGasType]);
    ];

cons = [
    boxCons;
    QdCons;
    nodalGasBalanceCons;
    nodalGasCompositionCons;
    pipeGasCompositionCons;
    ];
%
objfcn = objfcnUK(PGs,LCg,mpc) + penalty_Qd * sum(sum(slack_Qd))+ penalty_Qd * sum(sum(slack_nodal)) + 1*sum(sum(gasflow));

options = sdpsettings('verbose',2,'solver','ipopt', 'debug',1,'usex0',0);
options.gurobi.MIPgap = 1e-2; 
options.ipopt.tol = 1e-4; options.ipopt.max_iter = 1e6;
output = optimize(cons, objfcn, options);
%
% output = optimize(cons, objfcn);
%
% option1 = sdpsettings('verbose',2,'solver','gurobi', 'debug',1,'usex0',0,'relax',2);
% options.gurobi.MIPgap = 1e-2; options.ipopt.tol = 1e-4;
% output = optimize(cons, objfcn, option1);
% 
% gamma = value(gamma);

%%
PGs = value(PGs);
gasflow = value(gasflow);
gamma = value(gamma);
Qd = value(Qd);% Mm3/day
gasComposition = value(gasComposition); 
slack_Qd = value(slack_Qd);
LCg = value(LCg);
f_r = value(f_r);
objfcn = value(objfcn);
flowDirection = value(flowDirection);
slack_nodal = value(slack_nodal);
nodalgasInjection = value(nodalgasInjection);
gasCompositionCal = nodalgasInjection ./repmat(sum(nodalgasInjection,2),[1,nGasType]); 
save
% 可以试试先把gasflow的方向定下来
%% test
for r = 1:nGasType
    PGsbus = mpc.Gsou(:,1) ; 
    Cgs_PGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
    Qdbus = iGasLoad; 
    Cgs_Qd = sparse(Qdbus, (1:nGasLoad)', 1, nGb, nGasLoad); % connection matrix
    f_r(:,r) = Cgs_PGs*(PGs .* mpc.gasComposition(:,r)) - Cgs_Qd * Qd(:,r); % supply-demand, Mm3/day    
    % gas flow
    for  m = 1:nGl
        fb = mpc.Gline(m,1); tb = mpc.Gline(m,2);
        f_r(fb,r) = f_r(fb,r) - flowDirection(m) * gasflow(m,r);
        f_r(tb,r) = f_r(tb,r) + flowDirection(m) * gasflow(m,r);
    end
end
a=1
%% find shortest path
mpc.Gsou(10:11,2) = [151.88,15.9]/10;
for i = 1:nGb
    for n = 1:nGs
        gasSourceBus = mpc.Gsou(n,1);
        distance(i,n) = abs(gasSourceBus-i);
%         if distance(i,n) == 0
%             distance(i,n) = 1e-1;
%         end
    end
end
volume = mpc.Gsou(:,2)';
impact1 = repmat(volume,[nGb,1]) .* (160-distance);
weight = impact1./sum(impact1,2);

gasCompostionFake = zeros(nGb,nGasType);
for i = 1:nGb
    for r = 1:nGasType
        for n = 1:nGs
            gasCompostionFake(i,r) = gasCompostionFake(i,r) + weight(i,n)*mpc.gasComposition(n,r);
        end
    end
end