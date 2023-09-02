function W = nodalGasInjection(PGs,QptgInbus,gamma,gasFlow_sum,mpc,nGs,nGb,nGl)
PGsbus = mpc.Gsou(:,1) ; 
CgsPGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
W =  CgsPGs *PGs;% gas source term
W = W + QptgInbus; % PTG term
for i = 1:nGb
    for m = 1:nGl % get all the bus and Gline connected as injection
        if mpc.Gline(m,2) == i % i is the to bus
            % gamma = sign(gasFlow_sum),应该是这个关系
            W(i) = W(i) + (1+gamma(m))/2 * gasFlow_sum(m); 
        elseif mpc.Gline(m,1) == i % i is the from bus
            W(i) = W(i) + (-1+gamma(m))/2 * gasFlow_sum(m);
        end
    end
end

end