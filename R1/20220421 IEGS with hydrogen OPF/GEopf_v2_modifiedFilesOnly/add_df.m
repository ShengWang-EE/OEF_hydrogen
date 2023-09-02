function [add_df] = add_df(df,x,vv,mpc)
%ADD_DF renew df
%   Detailed explanation goes here
iPGs=vv.i1.PGs:vv.iN.PGs;
iLCg=vv.i1.LCg:vv.iN.LCg;
%%
gasCDF = calculateGasCDF(mpc);
df(iPGs) = mpc.Gcost;
df(iLCg) = gasCDF;
%
add_df=df;
end

