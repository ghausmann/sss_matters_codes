%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Auxiliary model (SSS recalibrated "rec" solution)
% Estimation of parameters (Table 1 in the paper)
% Note: This version also estimates the inverse of EIS, and sets a more
% conventional value for the persistence of TFP shocks
%---------------------------------------------------------------------------

clear;

%Choose whether to perform a full estimation (1), or just report moments (0)
perform_estimate = 0;
%Choose whether to save estimation results
save_estimate = 0;
%Choose whether to also compute Euler errors and the ergodic distribution of NFA.
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
model_path = genpath(fullfile(my_root,'model pac sss'));
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
etapar = 1001;
delta=0.014;            
alppha = 0.32;          
rho_x=0.99;             %CALIBRATION CHANGE
%rho_x=0.95;
rho_eps = 1;            %perturbation parameter persistence
%Auxiliary parameters
psibetta = 0;
psi_k = 0;
psi_i= 0; 

%Country-specific parameters (Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);
betta = 1/(1+(r_bar));   %discount factor

%To-be estimated parameters. Here just initialize them with some numbers
phi_d = 0.001; %controls Uzawa
phipar = 5; %controls capital adjustment costs
sigma_x = -3; %std. of innovations to productivity
D_bar = 2; %SSS debt
D1 = D_bar;
nu = 1;

%Make vector of parameter values and matrices of innovations
routine_rec_tasks;
eps_ind = find(x=='eps'); %index of epsilon in state vector xt

%--------------------------------------------------------------------------
% Estimation of the model
%--------------------------------------------------------------------------

%Generate matrix with shocks
Tt = 95;
Ti = 100;
T = Ti+Tt;
%Sigma = eye(5);
replications = 200;
%Vectors of param values to generate matrix with shocks
sP1 = [rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb];
sP2 = [rho_eps_r sigma_r_bar rho_sigma_r eta_r r_bar 0];
[innovations,~] = my_innovations(T,replications,sP1,sP2);

load('my_est_params_pac_rec'); %load estimated parameters
%set target for moments from the data, and initial guesses for estimated parameters
target = [4.77;1.31;3.81;1.78;3.47;-0.75];
if exist('my_est_params','var')
    guess_params = my_est_params;
else
    guess_params = [0.01 16.76 log(0.041) 2.56 1];
end

%Make monotonic change of the guess (ensures routine picks values within
%domain)
guess_params_opt = guess_params;
guess_params_opt(1) = log(guess_params_opt(1));
guess_params_opt(2) = (guess_params_opt(2))^0.5;
guess_params_opt(5) = (guess_params_opt(5))^0.5;

%Estimation function. It returns, for each targeted moment, the difference
%between target and simulation
myfun=@(x)my_smm_diff_function_sss_extrap(x,target,replications,innovations,model,params,my_dss_params,M,eps_ind,approx0,approx1,Ti,1);
options = optimoptions('lsqnonlin','Display','iter');

%First, compute simulated moments given the initial estimation guess, compare them
%with targets, and evaluate objective function (squared sum of differences)
sim_moments = target + myfun(guess_params_opt);
my_comp0 = [target sim_moments]
obj_fun = (sim_moments-target)'*(sim_moments-target)

%Then, execute the whole estimation
my_lb = [log(1e-5), 0, -Inf, -Inf, 1];
my_ub = [log(0.01), Inf, Inf, Inf, Inf];

if perform_estimate ==1
    t =clock;
    my_est = lsqnonlin(myfun,guess_params_opt,my_lb,my_ub,options)
    my_time1 = etime(clock,t)
else
    my_est = guess_params_opt;
end

%estimated parameters
phi_d = exp(my_est(1));  % Uzawa elasticity
phipar=(my_est(2))^2;       % capital adjustment costs
sigma_x = (my_est(3));  % std. of TFP shocks
D1 = (my_est(4));D_bar = D1; %SSS of debt
nu = (my_est(5))^2;       % inverse of EIS
% 
%Compute simulated moments from estimation, compare them
%with targets, and evaluate objective function (squared sum of differences)
sim_moments1 = target + myfun(my_est);
my_comp = [target sim_moments1]
obj_fun1 = (sim_moments1-target)'*(sim_moments1-target)
 
% %save the results
my_est_params = [phi_d phipar sigma_x D1 nu]

if save_estimate==1
    save_path = fullfile(my_root,'estimated parameters\my_est_params_pac_rec.mat');
    save(save_path,'my_est_params');
end

if compute_errors==1
    run fgru_pac_rec_errors.m;
end
