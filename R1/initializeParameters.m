function [GCVall, Mall, M_air, fsAll, aAll, R_air, T_stp, Prs_stp, Z, T_gas, eta_electrolysis,eta_methanation,etaGFU] = initializeParameters()
% CH4, C2H6, C3H8, C4H10, H2, N2, CO2
GCV_CH4 = 3.85 * 1e7;
GCV_C2H6 = 5.74 * 1e7;
GCV_C3H8 = 5.59 * 1e7;
GCV_C4H10 = 12.43 * 1e7;
GCV_hy = 12.75 * 1e6;      % J/m3
GCV_N2 = 0;     % J/m3
GCV_CO2 = 0;
GCVall = [GCV_CH4, GCV_C2H6, GCV_C3H8, GCV_C4H10, GCV_hy, GCV_N2, GCV_CO2];
GCV_ng = 41.04 * 1e6;     % J/m3
% 
M_CH4 = 16 * 1e-3;
M_C2H6 = 30 * 1e-3;
M_C3H8 = 44 * 1e-3;
M_C4H10 = 58 * 1e-3;
M_hy = 2 * 1e-3;           % kg/mol
M_N2 = 28 * 1e-3;
M_CO2 = 44 * 1e-3;
Mall = [M_CH4, M_C2H6, M_C3H8, M_C4H10, M_hy, M_N2, M_CO2];
M_ng = 17.478 * 1e-3;     % kg/mol
M_air = 29 * 1e-3;         % kg/mol
%  ﬂame speed 
fs_CH4 = 148;
fs_C2H6 = 301;
fs_C3H8 = 398;
fs_C4H10 = 513;
fs_hy = 339;           % kg/mol
fs_N2 = 0;
fs_CO2 = 0;
fsAll = [fs_CH4, fs_C2H6, fs_C3H8, fs_C4H10, fs_hy, fs_N2, fs_CO2];
% a
a_CH4 = 0.3;
a_C2H6 = 0.75;
a_C3H8 = 0.9759;
a_C4H10 = 1.0928;
a_hy = 0;         
a_N2 = 0.699;
a_CO2 = 0.9759;
aAll = [a_CH4, a_C2H6, a_C3H8, a_C4H10, a_hy, a_N2, a_CO2];
%
R_air = 287;               % J/(kg*K)
T_stp = 288;               % K
Prs_stp = 101325;          % Pa
Z = 1;                     % dimenssionless
T_gas = 281.15;            % K
eta_electrolysis = 0.7;    % from the energy perspective, the effciency is about 80%
eta_methanation = 0.8;
etaGFU = 0.4211;           % from the energy perspective, 从1/200换算而来
%%
WI_hy = GCV_hy / sqrt(M_hy/M_air);
WI_CH4 = GCV_CH4 / sqrt(M_CH4/M_air);
end