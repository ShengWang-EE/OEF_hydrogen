function [para] = initializeParameters2()
% p is the pressure (Pa), Z is the compressibility factor,
% Rspec is the specific gas constant for natural gas (J/kgK), 
% Θis the absolute temperature (K), ρ is the density (kg/m3),
% e the speed of wave propogation in the gas(m/s),
% M is the pipe flow rate (kg/s), 
% x is the distance along the pipeline (m), t is time (s), 
%d is the pipe diameter (m), A is the pipe’s crosssectional area (m2),
%F is the Fanning transmission factor(dimensionless)
% 气体的标准状况的气压温度是有标准规定的
% f=1/F^2, F is Fanning transmission factor, f is friction factor
% 可用稳态公式的Cij推出F: dp2/dx = Q^2/(Cij^2*L)
%%
para.Z = 0.8;
% para.R = 8.314;%J/mol*K
para.R = 8314.51/16;%Rgas
% para.Rgas = para.R/0.016 ;%J/(kg·K),M(CH4) = 16;
para.T = 281.15; %K
para.D = 0.89; % m
% para.D = 0.5; % m
para.A = (para.D/2)^2 * pi; %m2
para.L = 4000;%m pipeLength
para.e = sqrt(para.R*para.Z*para.T);

para.pn = 101325;% Pa, pressure in standard condition(STD) of gas
para.Tn = 298; %K, temerature in STD (in ref<multi...>is set as 288K
para.rhon = para.pn / (para.Z * para.R * para.Tn);
% Csquare = 9.07027/10^10*(10^6/24/3600)^2; % 注意单位换算p和Q的单位, for a representative pipe in ref......





end
