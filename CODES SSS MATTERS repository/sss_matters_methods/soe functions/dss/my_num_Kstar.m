function y = my_num_Kstar(my_dss_params,K1)

alppha = my_dss_params(1);
r_bar = my_dss_params(2);
delta = my_dss_params(3);
etapar = my_dss_params(4);
nu = my_dss_params(5);
D_bar = my_dss_params(6);
C_bar = my_dss_params(8);


Y_K = (r_bar+delta)/alppha;
H_K = Y_K^(1/(1-alppha));

Y = K1*Y_K;
H = K1*H_K;
C = C_bar*( (1-alppha)*Y/(H^(etapar)) )^(1/nu); 
I = (delta)*K1;

BC = Y - C - I + (D_bar/(1+r_bar) - D_bar);
%keyboard;
y = BC;

