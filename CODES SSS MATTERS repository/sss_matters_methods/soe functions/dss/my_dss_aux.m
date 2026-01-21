function [nxss,nyss,psi_k,C_bar] = my_dss_aux(sss_sol,my_dss_params)

D1 = sss_sol(1);
lK1 = sss_sol(2);
K1 = exp(lK1);
psi_i = sss_sol(3);

alppha = my_dss_params(1);
r_bar = my_dss_params(2);
delta = my_dss_params(3);
etapar = my_dss_params(4);
nu = my_dss_params(5);
sigma_r_bar = my_dss_params(6);
sigma_tb_bar = my_dss_params(7);

psi_k = alppha*( ( (1-alppha)^(1/(etapar-1)) )/( K1  ) )^( (etapar-1)*(1-alppha)/(alppha+etapar-1)   ) - (delta+(r_bar));
Y_K = (r_bar+delta+psi_k)/alppha;
H_K = Y_K^(1/(1-alppha));

Kstar = K1;
Ystar = K1*Y_K;
Hstar = K1*H_K;
Istar = (delta+psi_i)*K1;
Cstar = D1/(1+(r_bar))-D1+Ystar-Istar;
C_bar = Cstar;
lambdastar = (Cstar)^(-nu);
phistar = lambdastar;

%DSS for states
nxss=[D1;log(Kstar);log(Istar);sigma_r_bar;sigma_tb_bar;0;0;0;0;0;0];
%DSS for controls
nyss=[log(Cstar);log(Hstar);(r_bar);log(lambdastar);log(Ystar);log(phistar);log(Istar);D1];
