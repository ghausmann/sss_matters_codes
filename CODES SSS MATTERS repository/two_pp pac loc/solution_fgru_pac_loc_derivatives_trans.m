%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Location model (Appendix B), approximation around different points in the state
% space
% Computes numerical derivatives of log(consumption) w.r.t. an uncertainty shock
% Uses a THIRD order solution
%---------------------------------------------------------------------------

clear;

approx0 = 2; %order of approximation for SSS
approx1 = 3; %order of approximation for solution
% Choose algorithm for the Sylvester equation
algo='vectorize'; % Simple vectorization.
%algo='dlyap'; % Hessenberg-Schur algorithm.
%algo='gensylv'; % Kamenik algorithm

%Start fresh with the search path, and activate the required subfolders
my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
methods_path = genpath(fullfile(my_root,'sss_matters_methods'));
model_path = genpath(fullfile(my_root,'model pac loc'));
addpath(methods_path);
addpath(model_path);
addpath(fullfile(my_root,'estimated parameters'));

% Load the model:
load('model')

%Calibration 
%Common parameters:
rho_eps_tb=0.95;
sigma_tb_bar=-8.06;     %8.05 in paper, but 8.06 in code
rho_sigma_tb=0.94;
eta_tb=0.13;
nu=5;
etapar = 1001;
delta=0.014;            %depreciation
alppha = 0.32;          %capital income share
rho_x=0.95; 
rho_eps = 1;            %perturbation parameter persistence
%Auxiliary parameters
psi_k = 0;
psi_i= 0;
psi_d=4e-4;
Da = 0;

%Country-specific parameters (Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);

load('my_est_params_pac_dss_argentina'); %load estimated parameters
phi_d = my_est_params(1); %controls Uzawa
phipar = my_est_params(2); %controls capital adjustment costs
sigma_x = my_est_params(3); %std. of innovations to productivity
D1 = my_est_params(4); %SSS debt

D_bar = D1;
betta = 1/(1+(r_bar));   %discount factor

%Make vector of parameter values and matrices of innovations
routine_sss_loc_tasks;
eps_ind = find(x=='eps'); %index of epsilon in state vector xt

%Solve for the SSS
sss_sol = compute_sss_fgru_3_locp_num(model,params,my_dss_params,M,eta,eps_ind,approx0,[1 log(Kstar0) 0])

Dsss = sss_sol(1);
lK1 = sss_sol(2);
Ksss = exp(lK1);
psi_i = sss_sol(3);
Isss = (delta+psi_i)*Ksss;

%Make grids for the following variables: Debt, Capital, Investment, aux. parameter psi_d.
%The values in the grids imply a transition of the approximation point,
%from SSS to DSS. 
n_a = 25; %number of elements in the grids
Dmin = Dsss;
D_grid = [Dmin Dmin + ((1:(n_a-1))/(n_a-1)) * (D_bar - Dmin)];
Kmin = Ksss;
K_grid = [Kmin Kmin + ((1:(n_a-1))/(n_a-1)) * (Kstar0 - Kmin)];
Imin = Isss;
I_grid = [Imin Imin + ((1:(n_a-1))/(n_a-1)) * (Istar0 - Imin)];
psilsss = psi_d;
psildss = phi_d;
psil_grid = [psildss psildss + ((1:(n_a-1))/(n_a-1)) * (psilsss - psildss)];
psil_grid = psil_grid(end:-1:1);

my_h = 1e-5; %small increment of volatility for numerical derivative
%Preallocate vectors for derivatives
my_deriv = zeros(8,n_a);
deriv_cons = zeros(1,n_a);

for n=1:n_a
    
    %Set approximation point
    Da = D_grid(n);
    K1 = K_grid(n);
    Istar = I_grid(n);
    psi_i = (Istar/K1) - delta; 
    psi_d = psil_grid(n);

    [nxss,nyss,psi_k] = my_num_dss([Da log(K1) psi_i],my_dss_params);
    
    params(symparams=='psi_i') = psi_i;
    params(symparams=='psi_d') = psi_d;
    params(symparams=='Da') = Da;
    params(symparams=='psi_k') = psi_k;
    
    %Calculate derivatives of the policy functions abouts the DSS
    derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx1,algo);
    
    x0=nxss; % start at the steady state
    derivs.hx(eps_ind,eps_ind) = 1;
    x0(eps_ind) = 1; %impose the model of interest (epsilon=1)
    
    %Compute numerical derivatives
    x1 = x0; x1(x=='sigma_r')=x1(x=='sigma_r')+my_h;
    g0 = dr_gt(derivs,nyss,3,(x0-nxss));
    g1 = dr_gt(derivs,nyss,3,(x1-nxss));
    my_deriv(:,n) = (g1-g0)/my_h;
    deriv_cons(n) = my_deriv(1,n);

    
end

figure;plot(D_grid,deriv_cons);



