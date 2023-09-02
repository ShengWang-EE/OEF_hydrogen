function [h,WI] = consfcn_WobbeIndex(composition_hy,composition_gas,S_sqrt,GCV, nGb,GCV_gas,M_hy,M_gas,M_air)
WI0 = GCV_gas / sqrt(M_gas/M_air);
for i =1:nGb
    WI(i) = GCV(i) / S_sqrt(i);
    h(i) = WI(i) / WI0;
end
end