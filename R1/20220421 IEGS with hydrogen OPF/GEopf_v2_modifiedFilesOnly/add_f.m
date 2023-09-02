function [f_add_cost] = add_f(x,vv,mpc)
% calculate the additional cost
%   including gas purchasing cost, LCgcost (LCe cost, gas fired unit cost 
% are already included in the original cost function) 
%% 提取变量
PGs = x(vv.i1.PGs:vv.iN.PGs);
LCg = x(vv.i1.LCg:vv.iN.LCg);
%%
% noted the LCe cost are already calculated and updated in mpc
% the CDF of gas is setted as average of electricity CDF

% calculate gas CDF
gasCDF = calculateGasCDF(mpc);

%%
addCostPGs = PGs'*mpc.Gcost;
addCostLCg = sum(gasCDF .* LCg);
%% 
f_add_cost = addCostPGs + addCostLCg;%记得还要修改导数部分
end

