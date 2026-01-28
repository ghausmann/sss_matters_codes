function y = my_dss_psi_k(my_dss_params,D1,K1,psi_i,psi_k)

alppha = my_dss_params(1);
r_bar = my_dss_params(2);
delta = my_dss_params(3);
omega = my_dss_params(4);
nu = my_dss_params(5);
D_bar = my_dss_params(6);
Phi_d = my_dss_params(7);
Cstar = my_dss_params(8);


Y_K = (r_bar+delta+psi_k)/alppha;
%Y_H = Y_K^-(alppha/(1-alppha));
H_K = Y_K^(1/(1-alppha));

Y = K1*Y_K;
H = K1*H_K;
C = Cstar*( (1-alppha)*Y/(H^(omega)) )^(1/nu); 
I = (delta+psi_i)*K1;

BC = Y - C - I + (D1/(1+r_bar) - D1) - 0.5*Phi_d*(D1-D_bar)^2;
%keyboard;
y = BC;