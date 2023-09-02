function h = consfcn_gasFlow(gasFlow, mpc, nGl)
for i = 1:nGl
    gasFlowMax = mpc.Gline(i,5);
    h(i) = gasFlow(i) - gasFlowMax;

end