%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the DSS
% Computes simulated moments (to replicate numbers from Table 3 in the paper)
%---------------------------------------------------------------------------

clear;

%Choose whether to also compute Euler errors, deterministic simulation, and
%the ergodic distribution of NFA. (This will take a while for a fifth-order
%solution..).
compute_errors = 1;
%Choose order of approximation for solution (results will take a while for a fifth-order)
approx1 = 3; 
% Choose algorithm for the Sylvester equation
algo='vectorize'; % Simple vectorization.
%algo='dlyap'; % Hessenberg-Schur algorithm.
%algo='gensylv'; % Kamenik algorithm

%Start fresh with the search path, and activate the required subfolders
my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
methods_path = genpath(fullfile(my_root,'sss_matters_methods'));
model_path = genpath(fullfile(my_root,'model pac dss'));
addpath(methods_path);
addpath(model_path);
addpath(fullfile(my_root,'estimated parameters'));

rng(1) %fix seed
% Load the model:
load('model')

%Calibration 
%Common parameters:
rho_eps_tb=0.95;
sigma_tb_bar=-8.06;     
rho_sigma_tb=0.94;
eta_tb=0.13;
nu=5;                   
etapar = 1001;
delta=0.014;            
alppha = 0.32;          
rho_x=0.95;
rho_eps = 1;

%Country-specific parameters (for Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);

betta = 1/(1+(r_bar));   %discount factor

%To-be estimated parameters. Here just initialize them with some numbers
phi_d = 0.0006; %controls PAC
phipar=52; % capital adjustment costs
sigma_x=-2.9434; %std. innovations to productivity (in logs)
D_bar= 22.63; % DSS debt

%Perform routine DSS and moment-setting tasks
routine_est_dss_tasks;

%--------------------------------------------------------------------------
% Estimation of the model
%--------------------------------------------------------------------------

%Generate matrix with shocks
T0 = 3000; %initial number of periods in simulations (before innovations occur)
Tt = 95; %number of periods in each stochastic simulation
Ti = 100; %burn-in number of periods
T = Ti+Tt;
replications = 200;
sP1 = [rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb];
sP2 = [rho_eps_r sigma_r_bar rho_sigma_r eta_r r_bar 0];
[innovations,rates] = my_innovations(T,replications,sP1,sP2);

load('my_est_params_pac_dss_argentina');
if exist('my_est_params','var')
    guess_params = my_est_params;
else
    guess_params = [0.0001   64.8784   -3.2567   25.7968];
end

%Make monotonic change of the guess (ensures routine picks values within
%domain)
guess_params_opt = guess_params;
guess_params_opt(1) = log(guess_params_opt(1));
guess_params_opt(2) = (guess_params_opt(2))^0.5;

my_est = guess_params_opt;

%estimated parameters
phi_d = exp(my_est(1)); %controls PAC
phipar=my_est(2)^2; %capital adjustment costs
sigma_x = (my_est(3)); %std. innovations to productivity
D1 = (my_est(4)); %DSS of debt

%Report full set of moments of interest (would be nice to have a table
%here)
t =clock;
yf = my_smm_moments_dss_compare_pctile(my_est,replications,innovations,model,params,M,T0,Ti,approx1,0,95)
%yf = my_smm_moments_dss_compare_pctile(my_est,replications,innovations,model,params,M,T0,Ti,approx1,1,100)
my_time1 = etime(clock,t)
%ym = my_smm_moments_dss_compare_mod(my_est,replications,innovations,model,params,M)

my_est_params = [phi_d phipar sigma_x D1]

if compute_errors==1
    fgru_pac_dss_errors_pctile;
end