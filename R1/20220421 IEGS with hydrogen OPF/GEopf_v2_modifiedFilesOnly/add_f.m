function [f_add_cost] = add_f(x,vv,mpc)
% calculate the additional cost
%   including gas purchasing cost, LCgcost (LCe cost, gas fired unit cost 
% are already included in the original cost function) 
%% ��ȡ����
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
f_add_cost = addCostPGs + addCostLCg;%�ǵû�Ҫ�޸ĵ�������
end

