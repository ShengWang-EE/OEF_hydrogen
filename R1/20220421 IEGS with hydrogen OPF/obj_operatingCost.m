function [f,totalCost,genAndLCeCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy,penalty_PHI,penalty_sigma_PHI, penalty_sigma_x1000] = ...
    obj_operatingCost(Pg,PGs,LCg,Qptg,PHI, sigma_PHI, sigma_x1000, alpha_PHI,alpha_x, mpc,lambda)
% unit: $/hour
%

% 
CDFg = calculateGasCDF(mpc)*1;
Pg = mpc.baseMVA*Pg;
%

    genAndLCeCost = sum(Pg' * mpc.gencost(:,6));
%     sum(totcost_yalmip(mpc.gencost, Pg(:,:)')); % Pg includes GFU, but the gencost are 0.
    gasPurchasingCost = sum(PGs' * mpc.Gcost);
    gasCurtailmentCost = sum(LCg) * CDFg;

% subsidy of hydrogen and methane productions
% totalCost = totalCost- 0.1/6*1e6*sum(sum(Qptg));  % original
PTGsubsidy = sum(sum(Qptg))  *1e6/24 * 0.089 * 2.2/ 6.7;  
% penalty for gas flow
penalty_PHI = lambda * sum(PHI);
% penalty for gas flow concave
penalty_sigma_PHI = alpha_PHI * sum(sigma_PHI);
% penalty for gas composition
penalty_sigma_x1000 = alpha_x * sum(sum(sigma_x1000));

%%
totalCost = genAndLCeCost + gasPurchasingCost + gasCurtailmentCost - 1*PTGsubsidy;
f  = genAndLCeCost + gasPurchasingCost + gasCurtailmentCost - 1*PTGsubsidy + 0.1 * penalty_PHI + 1*penalty_sigma_PHI + 100*penalty_sigma_x1000;
end