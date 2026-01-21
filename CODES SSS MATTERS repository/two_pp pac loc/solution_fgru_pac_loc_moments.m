%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Auxiliary model (SSS fixed), approximation around the SSS
% Computes simulated moments (Tables 1, 3, and 5 in paper)
%---------------------------------------------------------------------------

clear;
%Choose if you want a subset of the replications (pctile of net exports)
%Default is pctile=100 (all replications)
pctile = 100;
%Choose if you want pruning active (default is yes)
act_pruning = 1;
%Choose if you want to activate a robustnetss check (default is not) 
check_phi_d = 0; %fixed PAC parameter
check_nu = 0; %log utility
%Choose approximation order
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

rng(1)
% Load the model:
load('model')

%Calibration 
%Common parameters:
rho_eps_tb=0.95;
sigma_tb_bar=-8.06;     
rho_sigma_tb=0.94;
eta_tb=0.13;
if check_nu==1
    nu=1;                   
else
    nu=5;
end
etapar = 1001;
delta=0.014;            
alppha = 0.32;          
rho_x=0.95;            
rho_eps = 1;            
%Auxiliary parameters
psi_k = 0;
psi_i= 0; 
Da = 0;
%Set auxiliary PAC parameter
if check_nu==1
    psi_d = 0.0014;
elseif check_phi_d==1
    psi_d=6e-4;
else
    psi_d=4e-4;
end

%Country-specific parameters (Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);

if check_nu==1
    load('my_est_params_pac_dss_argentina_nu');
elseif check_phi_d==1
    load('my_est_params_pac_dss_argentina_phi_d');
else
    load('my_est_params_pac_dss_argentina'); %load estimated parameters
end
    
phi_d = my_est_params(1); %controls Uzawa
phipar = my_est_params(2); %controls capital adjustment costs
sigma_x = my_est_params(3); %std. of innovations to productivity
D1 = my_est_params(4); %SSS debt

D_bar = D1;
betta = 1/(1+(r_bar));   %discount factor

%Make vector of parameter values and matrices of innovations
routine_sss_loc_tasks;
eps_ind = find(x=='eps'); %index of epsilon in state vector xt

%--------------------------------------------------------------------------
% Estimation of the model
%--------------------------------------------------------------------------

%Generate matrix with shocks
Tt = 95;
Ti = 100;
T = Ti+Tt;
replications = 200;
sP1 = [rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb];
sP2 = [rho_eps_r sigma_r_bar rho_sigma_r eta_r r_bar 0];
[innovations,~] = my_innovations(T,replications,sP1,sP2);

%Make monotonic change of the guess (ensures routine picks values within
%domain)
guess_params_opt = my_est_params;
guess_params_opt(1) = log(guess_params_opt(1));
guess_params_opt(2) = (guess_params_opt(2))^0.5;
my_est = guess_params_opt;

%Report full set of moments
t =clock;
yf = my_smm_moments_sss_compare_locp(my_est,replications,innovations,model,params,my_dss_params,M,eps_ind,approx0,approx1,Ti,act_pruning,pctile)
my_time1 = etime(clock,t)

fgru_pac_loc_errors;
