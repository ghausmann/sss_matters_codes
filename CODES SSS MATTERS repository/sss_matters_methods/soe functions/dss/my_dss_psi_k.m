function y = my_dss_psi_k(my_dss_params,D1,K1,psi_i,psi_k)

alppha = my_dss_params(1);
r_bar = my_dss_params(2);
delta = my_dss_params(3);
etapar = my_dss_params(4);
nu = my_dss_params(5);
D_bar = my_dss_params(6);
phi_d = my_dss_params(7);
C_bar = my_dss_params(8);


Y_K = (r_bar+delta+psi_k)/alppha;
H_K = Y_K^(1/(1-alppha));

Y = K1*Y_K;
H = K1*H_K;
C = C_bar*( (1-alppha)*Y/(H^(etapar)) )^(1/nu); 
I = (delta+psi_i)*K1;

BC = Y - C - I + (D1/(1+r_bar) - D1) - 0.5*phi_d*(D1-D_bar)^2;
%keyboard;
y = BC;

