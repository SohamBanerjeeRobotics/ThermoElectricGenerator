clc; clear; close all;
Power = zeros(1, 91); Efficiency = zeros(1, 91);
%% Parameters of the TEG Model. Bi2TE3 is taken as Semi-conductor.
syms I T_H T_C l_sc l_g n_NP;
% j = current density vector.
% T_H = Hot side temperature = 500 K
% T_C = Cold side temperature = 313 K
% l_sc = length of the semiconductor, same for both p and n type = 7.5 mm.
% l_g = distance between two legs = 1.5 mm
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
syms l_n k_n alpha_n mu_n rho_n
% l_n = width of the n-type semiconductor = 1.4 mm
% k_n = thermal conductivity of n-type semicoductor = 1.66 W/m.K
% alpha_n = Seebeck coefficient of n-type semiconductor = 190 * 10^(-6) V/K
% mu_n = Thompson coefficient of n-type semiconductor = 0.8 * 10^(-6) V/K
% rho_n = specific resistance of n-type semiconductor = 1.5 * 10^(-5) ohm m. 
A_n = ((-1) * (mu_n + alpha_n) * (I/l_n))/k_n; B_n = (rho_n * ((I/l_n) ^ 2))/k_n;
m1_n = ((-A_n) + sqrt((A_n^2) - (4 * B_n)))/2; 
m2_n = ((-A_n) - sqrt((A_n^2) - (4 * B_n)))/2;

% Parameters for the p-type semiconductor.
syms l_p k_p alpha_p mu_p rho_p
% l_p = width of the p-type semiconductor = 1.4 mm
% k_p = thermal conductivity of p-type semicoductor = 1.66 W/m.K
% alpha_p = Seebeck coefficient of p-type semiconductor = -190 * 10^(-6) V/K 
% mu_p = Thompson coefficient of p-type semiconductor = -0.8 * 10^(-6) V/K
% rho_p = specific resistance of p-type semiconductor = 1.5 * 10^(-5) ohm m. 
A_p = ((-1) * (mu_p + alpha_p) * (I/l_p))/k_p; B_p = (rho_p * ((I/l_p) ^ 2))/k_p; % modified
m1_p = ((-A_p) + sqrt((A_p^2) - (4 * B_p)))/2; 
m2_p = ((-A_p) - sqrt((A_p^2) - (4 * B_p)))/2;

% Parameters for the copper plate on the cold side.
% Same as hot side.

% Parameters for the ceramic plate on the cold side.
% Same as hot side.

% Parameters of external circuit.
syms R_L
% R_L = Load in the external circuit.

%% TEG Model
% We are going to solve an algebraic equation.
Coefficient_Matrix = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, ((-k_cu)/k_ce), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, (l_cu * (l_p + l_n)), (l_p + l_n), 0, (((((-k_p) * l_p)/k_cu) * m1_p) - ((I * alpha_p)/k_cu)), (((((-k_p) * l_p)/k_cu) * m2_p) - ((I * alpha_p)/k_cu)), (((((-k_n) * l_n)/k_cu) * m1_n) + ((I * alpha_n)/k_cu)), (((((-k_n) * l_n)/k_cu) * m2_n) + ((I * alpha_n)/k_cu)), 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, (m1_p * exp(m1_p * l_sc)), (m2_p * exp(m2_p * l_sc)), ((k_n * l_n)/(k_p * l_p) * m1_n * exp(m1_n * l_sc)), ((k_n * l_n)/(k_p * l_p) * m2_n * exp(m2_n * l_sc)), 0, ((-k_cu)/(k_p * l_p) * (l_p + l_n)), (-((I * alpha_n)/(k_p * l_p)) + ((I * alpha_p)/(k_p * l_p))), 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, l_cu, 1, 0, ((-k_ce)/k_cu), 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, l_ce, 1;
    l_ce, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, ((l_cu) ^ 2)/2, l_cu, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, ((l_cu) ^ 2)/2, l_cu, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, exp(m1_p * l_sc), exp(m2_p * l_sc), 0, 0, 0, 0, -1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, exp(m1_n * l_sc), exp(m2_n * l_sc), 0, 0, -1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, ((l_cu) ^ 2)/2, l_cu, 1, 0, -1;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
Output_Vector = [T_H; 0; 0; 0; 0; T_C; 0; 0; 0; 0; 0; 0; ((-rho_cu) * ((I/l_cu) ^ 2)) / k_cu; ((-rho_cu) * ((I/l_cu) ^ 2)) / k_cu];

% Variables : C1_ce1, C2_ce1, C1_cu1, C2_cu1, C3_cu1, C1_p, C2_p, C1_n,
% C2_n, C1_cu2, C2_cu2, C3_cu2, C1_ce2, C2_ce2.
Variable_vector = inv(Coefficient_Matrix) * Output_Vector; %#ok<MINV> 

% Getting the current density.
Var_eqn_LHS = ((-k_ce) * n_NP * (l_p + l_g + l_n) * Variable_vector(1,1)) + (k_ce * n_NP * (l_p + l_g + l_n) * Variable_vector(13,1)) - ((I^2) * R_L);

counter = 1;
for Load_R = 1 : 0.1 : 10 % varying the load resistance from 50 to 250 ohm.

eqn_LHS = subs(Var_eqn_LHS, [R_L, T_C, T_H, alpha_n, alpha_p, k_ce, k_cu, k_n, k_p, l_ce, l_cu, l_g, l_n, l_p, l_sc, mu_n, mu_p, n_NP, rho_cu, rho_n, rho_p], [Load_R, 296, 335, 207.135 * (10 ^ (-6)), (-207.135) * (10 ^ (-6)), 142, 142, 1.5205, 1.5205, 0.0091, 0.010, 0.0004, 0.0008, 0.0008, 0.0025, 0, 0, 125, 10.13, 10.645 * (10 ^ (-6)), 10.645 * (10 ^ (-6))]);
    
% Solving Equation using N-R method.
I_Initial = 20; % Taking the initial value of current to be 2 amps.
d_eqn_LHS = diff(eqn_LHS, I); % Differentiating the L.H.S of the equation.

for j = 1 : 25
    I_sol = I_Initial - (double(subs(eqn_LHS, I, I_Initial))/double(subs(d_eqn_LHS, I, I_Initial)));
    I_Initial = I_sol;
end
I_sol_real = real(I_sol); % Extracting the real part.
Power(1,counter) = (I_sol_real ^ 2) * Load_R; % Calculating the power.
Heat_In = real(double(subs(((-k_ce) * n_NP * (l_p + l_g + l_n) * Variable_vector(1,1)), [I, T_C, T_H, alpha_n, alpha_p, k_ce, k_cu, k_n, k_p, l_ce, l_cu, l_g, l_n, l_p, l_sc, mu_n, mu_p, n_NP, rho_cu, rho_n, rho_p],[I_sol_real, 296, 335, 207.135 * (10 ^ (-6)), (-207.135) * (10 ^ (-6)), 142, 142, 1.5205, 1.5205, 0.0091, 0.010, 0.0004, 0.0008, 0.0008, 0.0025, 0, 0, 125, 10.13, 10.645 * (10 ^ (-6)), 10.645 * (10 ^ (-6))])));
Efficiency(1,counter) = (Power(1,counter) / (Heat_In)) * 100;
counter = counter + 1;
end
%% Plotting
subplot(2,1,1);
plot(1 : 0.1 : 10, Power, '-b');
xlabel("Load Resistance (Ohm) --->");
ylabel("Output Power (W) --->");
title("Output Power vs Load Resistance");
grid on;

subplot(2,1,2);
plot(1 : 0.1 : 10, Efficiency, '-r');
xlabel("Load Resistance (Ohm) --->");
ylabel("Efficiency (%) --->");
title("Efficiency vs Load Resistance");
grid on;