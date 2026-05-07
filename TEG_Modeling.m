clc; clear; close all;
Power = zeros(1, 51); Efficiency = zeros(1, 51);
%% Parameters of the TEG Model. Bi2TE3 is taken as Semi-conductor.
syms I T_H T_C l_sc n_NP;
% j = current density vector.
% T_H = Hot side temperature.
% T_C = Cold side temperature = 313 K
% l_sc = length of the semiconductor, same for both p and n type = 7.5 mm.
% l_g = distance between two legs = 0.05 approx 0 mm
% n_NP = number of NP pairs = 150

% Parameters for the ceramic plate on the hot side.
syms l_ce k_ce; 
% l_ce = thickness of the ceramic plate = 0.5 mm 
% k_ce = thermal conductivity of ceramic = 230 W/m.K

% Parameters for the copper plate on the hot side.
syms l_cu k_cu rho_cu; 
% l_cu = thickness of copper plate = 3 mm 
% k_cu = thermal conductivity of copper plate = 398 W/m.K
% rho_cu = specific resistance of copper = 1.68 * 10^(-8) ohm m.

% Parameters for the n-type semiconductor.
syms l_n k_n alpha_n rho_n
% l_n = width of the n-type semiconductor = 1.4 mm
% k_n = thermal conductivity of n-type semicoductor = 1.66 W/m.K
% alpha_n = Seebeck coefficient of n-type semiconductor = -190 * 10^(-6) V/K
% * mu_n = Thompson coefficient of n-type semiconductor = 0.8 * 10^(-6) V/K
% rho_n = specific resistance of n-type semiconductor = 1.5 * 10^(-5) ohm m. 
A_n = ((-1) * alpha_n * (I/l_n))/k_n; B_n = (rho_n * ((I/l_n) ^ 2))/k_n;

% Parameters for the p-type semiconductor.
syms l_p k_p alpha_p rho_p
% l_p = width of the p-type semiconductor = 1.4 mm
% k_p = thermal conductivity of p-type semicoductor = 1.66 W/m.K
% alpha_p = Seebeck coefficient of p-type semiconductor = 190 * 10^(-6) V/K 
% * mu_p = Thompson coefficient of p-type semiconductor = -0.8 * 10^(-6) V/K
% rho_p = specific resistance of p-type semiconductor = 1.5 * 10^(-5) ohm m. 
A_p = ((-1) * alpha_p * (I/l_p))/k_p; B_p = (rho_p * ((I/l_p) ^ 2))/k_p; % modified

% Parameters for the copper plate on the cold side.
% Same as hot side.

% Parameters for the ceramic plate on the cold side.
% Same as hot side.

% Parameters of external circuit.
syms R_L
% R_L = Load in the external circuit = 150 Ohm

%% TEG Model
% We are going to solve an algebraic equation.
Coefficient_Matrix = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, ((-k_cu)/k_ce), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, (-k_cu) * l_cu * (l_p + l_n), (-k_cu) * (l_p + l_n), 0, ((k_p * l_p) + ((I * alpha_p)/A_p)), ((-1) * I * alpha_p), ((k_n * l_n) + (((-1) * I * alpha_n)/A_n)), (I * alpha_n), 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, (((k_p * l_p) + ((I * alpha_p)/A_p)) * exp((-1) * A_p * l_sc)), ((-1) * I * alpha_p), (((k_n * l_n) - ((I * alpha_n)/A_n)) * exp((-1) * A_n * l_sc)), (I * alpha_n), 0, ((-k_cu) * (l_p + l_n)), 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, (l_cu), 1, 0, (-1) * (k_ce/k_cu), 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, l_ce, 1;
    l_ce, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, ((l_cu) ^ 2)/2, l_cu, 1, (1/A_p), (-1), 0, 0, 0, 0, 0, 0, 0;
    0, 0, ((l_cu) ^ 2)/2, l_cu, 1, 0, 0, (1/A_n), (-1), 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, (exp((-1) * A_p * l_sc))/A_p, (-1), 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, (exp((-1) * A_n * l_sc))/A_n, (-1), 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, ((l_cu) ^ 2)/2, l_cu, 1, 0, (-1);
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
Output_Vector = [T_H; 0; (((k_p * l_p * B_p)/A_p) + ((k_n * l_n * B_n)/A_n)); (((k_p * l_p * B_p)/A_p) + ((k_n * l_n * B_n)/A_n) - ((alpha_p * I * B_p * l_sc)/A_p) + ((alpha_n * I * B_n * l_sc)/A_n)); 0; T_C; 0; 0; 0; ((-1) * B_p * l_sc)/A_p; ((-1) * B_n * l_sc)/A_n; 0; ((-rho_cu) * ((I/l_cu) ^ 2)) / k_cu; ((-rho_cu) * ((I/l_cu) ^ 2)) / k_cu];

% Variables : C1_ce1, C2_ce1, C1_cu1, C2_cu1, C3_cu1, C1_p, C2_p, C1_n,
% C2_n, C1_cu2, C2_cu2, C3_cu2, C1_ce2, C2_ce2.
Variable_vector = inv(Coefficient_Matrix) * Output_Vector; %#ok<MINV> 

% Getting the current density.
Var_eqn_LHS = ((-k_ce) * n_NP * (l_p + l_n) * Variable_vector(1,1)) + (k_ce * n_NP * (l_p + l_n) * Variable_vector(13,1)) + (alpha_n * I * ((((-1) * Variable_vector(8,1))/A_n) + Variable_vector(9,1))) - (alpha_p * I * ((((-1) * Variable_vector(6,1))/A_p) + Variable_vector(7,1))) - ((I^2) * R_L);
Var_eqn_LHS = Var_eqn_LHS + (alpha_p * I * (Variable_vector(7,1) - ((B_p * l_sc)/A_p) - ((Variable_vector(6,1) * exp((-A_p) * l_sc))/A_p))) - (alpha_n * I * (Variable_vector(9,1) - ((B_n * l_sc)/A_n) - ((Variable_vector(8,1) * exp((-A_n) * l_sc))/A_n)));

counter = 1;
for Hot_Side_Temp = 400 : 8 : 800 % varying the load resistance from 50 to 250 ohm.

eqn_LHS = subs(Var_eqn_LHS, [R_L, T_C, T_H, alpha_n, alpha_p, k_ce, k_cu, k_n, k_p, l_ce, l_cu, l_n, l_p, l_sc, n_NP, rho_cu, rho_n, rho_p], [150, 313, Hot_Side_Temp, (-190) * (10 ^ (-6)), 190 * (10 ^ (-6)), 230, 398, 1.66, 1.66, 0.0005, 0.003, 0.0014, 0.0014, 0.0075, 150, 1.68 * (10 ^ (-8)), 1.5 * (10 ^ (-8)), 1.5 * (10 ^ (-8))]);
    
% Solving Equation using N-R method.
I_Initial = 20; % Taking the initial value of current to be 2 amps.
d_eqn_LHS = diff(eqn_LHS, I); % Differentiating the L.H.S of the equation.

for j = 1 : 25
    I_sol = I_Initial - (double(subs(eqn_LHS, I, I_Initial))/double(subs(d_eqn_LHS, I, I_Initial)));
    I_Initial = I_sol;
end
I_sol_real = real(I_sol); % Extracting the real part.
Power(1,counter) = (I_sol_real ^ 2) * 150; % Calculating the power.
Heat_In = real(double(subs(((-k_ce) * n_NP * (l_p + l_n) * Variable_vector(1,1)), [I, T_C, T_H, alpha_n, alpha_p, k_ce, k_cu, k_n, k_p, l_ce, l_cu, l_n, l_p, l_sc, n_NP, rho_cu, rho_n, rho_p],[I_sol_real, 313, Hot_Side_Temp, (-190) * (10 ^ (-6)), 190 * (10 ^ (-6)), 230, 398, 1.66, 1.66, 0.0005, 0.003, 0.0014, 0.0014, 0.0075, 150, 1.68 * (10 ^ (-8)), 1.5 * (10 ^ (-8)), 1.5 * (10 ^ (-8))])));
Efficiency(1,counter) = (Power(1,counter) / (Heat_In)) * 100;
counter = counter + 1;
end
%% Plotting
subplot(2,1,1);
plot(400 : 8 : 800, Power, '-b');
xlabel("Hot side temperature (K) --->");
ylabel("Output Power (W) --->");
title("Output Power vs Hot side temperature");
grid on;

subplot(2,1,2);
plot(400 : 8 : 800, Efficiency, '-r');
xlabel("Hot side temperature (K) --->");
ylabel("Efficiency (%) --->");
title("Efficiency vs Hot side temperature");
grid on;