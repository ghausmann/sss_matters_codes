function [nxss,nyss] = my_num_dss_K(my_dss_params)
%This function computes the DSS for a perturbation model used as a initial
%guess for Taylor projection.

alppha = my_dss_params(1);
r_bar = my_dss_params(2);
delta = my_dss_params(3);
etapar = my_dss_params(4);
nu = my_dss_params(5);
D_bar = my_dss_params(6);
% phi_d = my_dss_params(7);
C_bar = my_dss_params(8);
sigma_r_bar = my_dss_params(9);
sigma_tb_bar = my_dss_params(10);


my_dss_fun = @(x)my_num_Kstar(my_dss_params,x);
options = optimoptions('fsolve','Display','off'); 
K_guess = my_Kstar(my_dss_params);
K1 = fsolve(my_dss_fun,K_guess,options);

Y_K = (r_bar+delta)/alppha;
H_K = Y_K^(1/(1-alppha));

Kstar = K1;
Ystar = K1*Y_K;
Hstar = K1*H_K;
Cstar = C_bar*( (1-alppha)*Ystar/(Hstar^(etapar)) )^(1/nu);
Istar = (delta)*K1;
lambdastar = (Cstar)^(-nu);
phistar = lambdastar;

%DSS for states
nxss=[D_bar;log(Kstar);log(Istar);sigma_r_bar;sigma_tb_bar;0;0;0];
%DSS for controls
nyss=[log(Cstar);log(Hstar);(r_bar);log(phistar);log(Istar);D_bar];
