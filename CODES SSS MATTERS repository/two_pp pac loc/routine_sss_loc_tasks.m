%This script performs routine tasks of making
%vectors of parameter values and matrices for innovations

%DSS of ORIGINAL model (as implied by DSS calibration)
Hstar0 = ( (1-alppha)*(alppha/((r_bar)+(delta)))^(alppha/(1-alppha))  )^(1/(etapar-1));
Kstar0 = Hstar0*(alppha/((r_bar)+(delta)))^(1/(1-alppha));
Istar0 = (delta)*Kstar0;
Ystar0 = (Kstar0^alppha)*Hstar0^(1-alppha);
Cstar0 = D_bar/(1+(r_bar))-D_bar+Ystar0-Istar0;
C_bar = Cstar0;
my_dss_params = [alppha r_bar delta etapar nu D_bar phi_d C_bar sigma_r_bar sigma_tb_bar];

% Make a vector of parameter values:
params=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, rho_eps, etapar, phi_d, C_bar, psi_d, psi_i, psi_k, Da];
% Make eta
% Define eta as in Schmitt-Grohe and Uribe (2004):
eta=[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;eta_r 0 0 0 0;0 eta_tb 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 exp(sigma_x);0 0 0 0 0];

% Compute the cross moments:
n_e=size(eta,2); % number of shocks.
M2f = eye(n_e);
M2f(1,1) = 0.00001; %This makes uncertainty shocks unexpected (only relevant for fifth-order solutions)
M2f(2,2) = 0.00001;
[n_nodes,nodes,weights] = Monomials_2(n_e,M2f); % this quadrature function was written by Judd, Maliar, Maliar and Valero (2014).
nodes=nodes'; % transpose to n_e-by-n_nodes
M=get_moments(nodes,weights,approx1);
