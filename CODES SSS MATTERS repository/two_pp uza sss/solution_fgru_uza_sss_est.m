%---------------------------------------------------------------------------
% FGRU Small open-economy model
% UZAWA model, approximation around the SSS
% Estimation of parameters (Table 5 in the paper)
%---------------------------------------------------------------------------

clear;

%Choose if you want pruning active (default is yes)
act_pruning = 1;
%Choose whether to perform a full estimation (1), or just report moments (0)
perform_estimate = 0;
%Choose whether to save estimation results
save_estimate = 0;
%Choose whether to also compute Euler errors
compute_errors = 1;
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
model_path = genpath(fullfile(my_root,'model uza sss'));
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
%Auxiliary parameters
psibetta = 0;
psi_k = 0;
psi_i= 0; 

%Country-specific parameters (Argentina)
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);

%To-be estimated parameters. Here just initialize them with some numbers
phi_d = 0.0782;
phipar=34.8722; 
sigma_x=-2.9434;
D_bar = 2; %SSS debt
D1 = D_bar;
betta = 1/(1+(r_bar));   %discount factor

%Make vector of parameter values and matrices of innovations
routine_sss_tasks;
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
 
load('my_est_params_uza_sss');
%set target for moments from the data, and initial guesses for estimated parameters
target = [(4.77);(1.31);(3.81);1.78]; %Target moments from Argentina data
if exist('my_est_params','var')
    guess_params = my_est_params;
else
    guess_params = [0.01 16.76 log(0.041) 2.56];
end

%Make monotonic change of the guess (ensures routine picks values within
%domain)
guess_params_opt = guess_params;
guess_params_opt(1) = log(guess_params_opt(1));
guess_params_opt(2) = (guess_params_opt(2))^0.5;

%Estimation function. It returns, for each targeted moment, the difference
%between target and simulation
myfun=@(x)my_smm_diff_function_sssp(x,target,replications,innovations,model,params,my_dss_params,M,eps_ind,approx0,approx1,Ti,act_pruning);
options = optimoptions('lsqnonlin','Display','iter');

%First, compute simulated moments given the initial estimation guess, compare them
%with targets, and evaluate objective function (squared sum of differences)
t =clock; 
sim_moments = target + myfun(guess_params_opt);
my_comp = [target sim_moments]
obj_fun = (sim_moments-target)'*(sim_moments-target)
my_time1 = etime(clock,t)

%Then, execute the whole estimation
%Lower and upper bounds for parameters
my_lb = [log(1e-05), 0, -Inf, -Inf];
my_ub = [log(0.01), Inf, Inf, Inf];

if perform_estimate ==1
    t =clock;
    my_est = lsqnonlin(myfun,guess_params_opt,my_lb,my_ub,options)
    my_time1 = etime(clock,t)
else
    my_est = guess_params_opt;
end

%estimated parameters
phi_d = exp(my_est(1)); %controls PAC
phipar=my_est(2)^2; %capital adjustment costs
sigma_x = (my_est(3)); %std. innovations to productivity
D1 = (my_est(4));D_bar = D1; %SSS of debt

yf = my_smm_moments_sss_comparep(my_est,replications,innovations,model,params,my_dss_params,M,eps_ind,approx0,approx1,Ti,act_pruning)

%save the results
my_est_params = [phi_d phipar sigma_x D1]

if save_estimate==1
    save_path = fullfile(my_root,'estimated parameters\my_est_params_uza_sss.mat');
    save(save_path,'my_est_params');
end

if compute_errors==1
    fgru_uza_sss_errors;
end