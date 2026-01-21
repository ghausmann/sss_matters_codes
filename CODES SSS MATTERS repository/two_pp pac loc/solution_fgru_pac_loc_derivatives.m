%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Auxiliary model (SSS fixed), approximation around the SSS
% Computes consumption elasticities about the SSS (Table 4 in the paper)
%---------------------------------------------------------------------------

clear;

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
approx1 = 5;
routine_sss_loc_tasks;
eps_ind = find(x=='eps'); %index of epsilon in state vector xt

%--------------------------------------------------------------------------
% Derivatives 3rd order solution at SSS
%--------------------------------------------------------------------------

%Solve for the SSS
sss_sol_3 = compute_sss_fgru_3_locp_num(model,params,my_dss_params,M,eta,eps_ind,2,[1 log(Kstar0) 0])
%Recompute DSS of auxiliary model
[nxss,nyss,psi_k3] = my_num_dss(sss_sol_3,my_dss_params);

params(symparams=='Da') = sss_sol_3(1);
params(symparams=='psi_k') = psi_k3;
params(symparams=='psi_i') = sss_sol_3(3);

%Calculate derivatives of the policy functions abouts the DSS
derivs=solve_dsge(model,params,M,eta,nxss,nyss,3,algo);

x0=nxss; % start at the steady state
derivs.hx(eps_ind,eps_ind) = 1;
x0(eps_ind) = 1; %impose the model of interest (epsilon=1)

my_h = 1e-5;
x1 = x0; x1(x=='sigma_r')=x1(x=='sigma_r')+my_h;
g0_3 = dr_gt(derivs,nyss,3,(x0-nxss));
g1_3 = dr_gt(derivs,nyss,3,(x1-nxss));
my_deriv_3 = (g1_3-g0_3)/my_h

%--------------------------------------------------------------------------
% Derivatives 5th order solution at SSS
%--------------------------------------------------------------------------

%Solve for the SSS
sss_sol_4 = compute_sss_fgru_3_locp_num(model,params,my_dss_params,M,eta,eps_ind,4,[1 log(Kstar0) 0])
%Recompute DSS of auxiliary model
[nxss,nyss,psi_k4] = my_num_dss(sss_sol_4,my_dss_params);

params(symparams=='Da') = sss_sol_4(1);
params(symparams=='psi_k') = psi_k4;
params(symparams=='psi_i') = sss_sol_4(3);

%Calculate derivatives of the policy functions abouts the DSS
derivs=solve_dsge(model,params,M,eta,nxss,nyss,5,algo);

x0=nxss; % start at the steady state
derivs.hx(eps_ind,eps_ind) = 1;
x0(eps_ind) = 1; %impose the model of interest (epsilon=1)

my_h = 1e-5;
x1 = x0; x1(x=='sigma_r')=x1(x=='sigma_r')+my_h;
g0_5 = dr_gt(derivs,nyss,5,(x0-nxss));
g1_5 = dr_gt(derivs,nyss,5,(x1-nxss));
my_deriv_5 = (g1_5-g0_5)/my_h;

%First row are consumption elasticities
compare_derivs_2param_loc = [my_deriv_3 my_deriv_5]
