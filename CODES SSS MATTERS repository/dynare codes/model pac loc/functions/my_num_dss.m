function [yss,psi_k] = my_num_dss(sss_sol,my_dss_params)

Da = sss_sol(1);
lK1 = sss_sol(2);
K1 = exp(lK1);
psi_i = sss_sol(3);

alppha = my_dss_params(1);
r_bar = my_dss_params(2);
delta = my_dss_params(3);
omega = my_dss_params(4);
nu = my_dss_params(5);
% D_bar = my_dss_params(6);
% psin = my_dss_params(7);
Cstar0 = my_dss_params(8);
sigma_r_bar = my_dss_params(9);
sigma_tb_bar = my_dss_params(10);


my_dss_fun = @(x)my_dss_psi_k(my_dss_params,Da,K1,psi_i,x);
options = optimoptions('fsolve','Display','off'); 
psi_k = fsolve(my_dss_fun,0,options);

Y_K = (r_bar+delta+psi_k)/alppha;
H_K = Y_K^(1/(1-alppha));

Kstar = K1;
Ystar = K1*Y_K;
Hstar = K1*H_K;
Cstar = Cstar0*( (1-alppha)*Ystar/(Hstar^(omega)) )^(1/nu);
Istar = (delta+psi_i)*K1;
lambdastar = (Cstar)^(-nu);
phistar = lambdastar;

%DSS for states
yss = [log(Cstar);log(Hstar);(r_bar);Da;...
    sigma_r_bar;sigma_tb_bar;0;0;0;Da;log(Kstar);log(Istar);0;...
    log(lambdastar);log(Ystar);log(phistar);log(Istar)];