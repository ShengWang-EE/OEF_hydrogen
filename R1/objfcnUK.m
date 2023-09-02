function objfcn = objfcnUK(PGs,LCg,mpc)
    CDF = 1e3;
    objfcn = sum(PGs(1:9))*mpc.Gcost + sum(LCg) * CDF;
end