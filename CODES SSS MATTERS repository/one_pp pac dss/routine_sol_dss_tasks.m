%This script performs the routine tasks of calculating the DSS, making
%vectors of parameter values and matrices for innovations, and computing the
%perturbation solution

%DSS of the model
Hstar = ( (1-alppha)*(alppha/((r_bar)+(delta)))^(alppha/(1-alppha))  )^(1/(etapar-1));
Kstar = Hstar*(alppha/((r_bar)+(delta)))^(1/(1-alppha));
Istar = delta*Kstar;
Ystar = (Kstar^alppha)*Hstar^(1-alppha);
Cstar = D1/(1+(r_bar))-D1+Ystar-Istar;
C_bar = Cstar;
lambdastar = (C_bar)^(-nu);
phistar = lambdastar;

nxss=[D1;log(Kstar);log(Istar);sigma_r_bar;sigma_tb_bar;0;0;0;0;0;0]; %DSS states
nyss=[log(Cstar);log(Hstar);(r_bar);log(lambdastar);log(Ystar);log(phistar);log(Istar);D1]; %DSS controls

% Make a vector of parameter values:
params=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, rho_eps, etapar, phi_d, C_bar];

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

% Compute the perturbation solution:
%DSS check (the DSS should make all equilibrium conditions zero)
mycheck0 = double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
mycheck=sum(abs(mycheck0))

%Calculate derivatives of the policy functions about the DSS
derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx1,algo);
