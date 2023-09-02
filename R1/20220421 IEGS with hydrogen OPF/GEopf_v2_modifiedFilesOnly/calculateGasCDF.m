function gasCDF = calculateGasCDF(mpc)
% calculate gas CDF 
dispatchableLoadIndex = mpc.originalGenNumber+1;
electricityCDF = mean(mpc.gencost(dispatchableLoadIndex,6));
Hg = 39;%MW/(m/s3) ���������ĵ�λ��Mm3/day
etaE = 1; etaG = 1;
gasCDF = electricityCDF * Hg * etaE / etaG * 1000000 /24 /3600;
% test---------------------------
% ����̬��ʱ��Ҳ���е�

end