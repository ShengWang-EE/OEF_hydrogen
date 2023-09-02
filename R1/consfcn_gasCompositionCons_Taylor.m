function g = consfcn_gasCompositionCons_Taylor(PGs0,Qptg0, gasFlow, composition_hy,composition_gas, mpc, signGf, nGb,nGl)
%% convert
PGs = PGs0 * 1e6/24/3600;
Qptg = Qptg0 * 1e6/24/3600;
g = [];

for i = 1:nGb
        connectedGlineIndex = []; connectedBus  =[];
        for m = 1:nGl % get all the bus and Gline connected as injection
            if (mpc.Gline(m,2) == i) 
                connectedGlineIndex = [connectedGlineIndex,m];
                connectedBus = [connectedBus, mpc.Gline(m,1)];
            elseif (mpc.Gline(m,1) == 1) 
                connectedGlineIndex = [connectedGlineIndex, m];
                connectedBus = [connectedBus, mpc.Gline(m,2)];
            end
        end

        ptgIndex = find(mpc.ptg(:,1)==i);

        if isempty(ptgIndex) == 1 % no ptg
            ptgInjection = 0;
        else
            ptgInjection = Qptg(ptgIndex);
        end

        flowInjection = 0; hydrogenFlowInjection = 0; naturalGasFlowInjection = 0;
        if isempty(fb) == 0 % has from bus
            for j = 1:size(fb,2)
                flowInjection = flowInjection + signGf(fromGlineIndex(j)) * gasFlow(fromGlineIndex(j));
                hydrogenFlowInjection = hydrogenFlowInjection + composition_hy(fb(j)) * signGf(fromGlineIndex(j)) * gasFlow(fromGlineIndex(j));       
                naturalGasFlowInjection = naturalGasFlowInjection + composition_gas(fb(j)) * signGf(fromGlineIndex(j)) * gasFlow(fromGlineIndex(j));   
            end
        else 
            flowInjection = 0;
            hydrogenFlowInjection = 0;
            naturalGasFlowInjection = 0;
        end

        if ismember(i,mpc.Gsou(:,1)) % has gas source
            gasSourceInjection = PGs(find(mpc.Gsou(:,1)==i));
        else % no gas source
            gasSourceInjection = 0;
        end

        g = [g;
            composition_hy(i) == (ptgInjection + hydrogenFlowInjection) / (ptgInjection + gasSourceInjection + flowInjection);
            composition_gas(i) == (gasSourceInjection + naturalGasFlowInjection) / (ptgInjection + gasSourceInjection + flowInjection);
            ];
end


end