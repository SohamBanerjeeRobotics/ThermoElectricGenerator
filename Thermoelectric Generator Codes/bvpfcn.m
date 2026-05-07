function dTdy = bvpfcn(~,T,rho_sc,K_sc,l,alpha_sc,alpha_sc_slope)
dTdy = zeros(2,1);
dTdy(1,1) = T(2,1);
dTdy(2,1) = (((-1) * rho_sc * I * I) / (K_sc * l * l)) + (T(1,1) * (alpha_sc_slope / K_sc) * (I / l) * T(2,1)) + ((alpha_sc / K_sc) * (I / l) * T(2,1));
end