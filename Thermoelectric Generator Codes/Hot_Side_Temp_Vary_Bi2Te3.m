%% 
% Initial assumptions for the variables : C1_Ce1, C1_Ce2, I

% Properties of Ceramic Substrate.

% Properties of Copper Plate.

% Properties of Semiconductor.

% From Top side.
% Ceramic 1: C1_Ce1, C2_Ce1
C2_Ce1 = T_H; % T_H = Hot side temperature.

% Copper 1: C1_Cu1, C2_Cu1, C3_Cu1
C1_Cu1 = ((-1) * rho_cu * (I ^ 2))/(K_cu * l_cu * l_cu);
C2_Cu1 = (K_ce / K_cu) * C1_Ce1; C3_Cu1 = (C1_Ce1 * l_ce) + C2_Ce1;

% From bottom side.
% Ceramic 2: C1_Ce2, C2_Ce2
C2_Ce2 = T_C - (C1_Ce2 * l_ce); % T_C = Cold side Temperature.

% Copper 2: C1_Cu2, C2_Cu2, C3_Cu2
C1_Cu2 = ((-1) * rho_cu * (I ^ 2))/(K_cu * l_cu * l_cu);
C2_Cu2 = ((K_ce / K_cu) * C1_Ce2) - (C1_Cu2 * l_cu); 
C3_Cu2 = C2_Ce2 - (C1_Cu2 * l_cu * l_cu * 0.5) - (C2_Cu2 * l_cu);

%% Solving the Boundary Value Problem.
Initial_T_SCP = (C1_Cu1 * (0.5 * l_cu * l_cu)) + (C2_Cu1 * l_cu) + C3_Cu1;
Initial_T_SCN = Initial_T_SCP;

Final_T_SCP = C3_Cu2; Final_T_SCN = Final_T_SCP;

solinit = bvpinit(linspace(0, l_sc, 100),[1,1]);
sol_p = bvp4c(@(y,T) bvpfcn(y,T,rho_p,K_p,l_p,alpha_p,alpha_p_slope), @(Ti,Tf) bcfcn(Ti,Tf,Initial_T_SCP,Final_T_SCP), solinit);
sol_n = bvp4c(@(y,T) bvpfcn(y,T,rho_n,K_n,l_n,alpha_n,alpha_n_slope), @(Ti,Tf) bcfcn(Ti,Tf,Initial_T_SCN,Final_T_SCN), solinit);
