%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the DSS
% Estimation of parameters
%---------------------------------------------------------------------------

clear;

%set country; you must not set more than one entry to 1
Argentina = 1;
Brazil = 0;
Ecuador = 0;
%Choose whether to perform a full estimation (1), or just report moments (0)
perfom_estimate = 0;
%Choose whether to save estimation results
save_estimate = 0;
%Choose whether to also compute Euler errors, deterministic simulation, and
%the ergodic distribution of NFA.
compute_errors = 1;
%Choose order of approximation for solution
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

%Country-specific parameters:
if Argentina ==1
    rho_eps_r=0.97;
    sigma_r_bar=-5.71;
    rho_sigma_r=0.94;
    eta_r=0.46;
    r_bar=(0.02);
    elseif Brazil ==1
    rho_eps_r=0.95;
    sigma_r_bar=-6.97;
    rho_sigma_r=0.95;
    eta_r=0.28;
    r_bar=(0.007);
    elseif Ecuador == 1
    rho_eps_r=0.95;
    sigma_r_bar=-6.06;
    rho_sigma_r=0.96; 
    eta_r=0.35;
    r_bar=0.011;
      
end

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
[innovations,~] = my_innovations(T,replications,sP1,sP2);

if Argentina == 1
    load('my_est_params_pac_dss_argentina');
elseif Brazil == 1
    load('my_est_params_pac_dss_brazil');
elseif Ecuador == 1
    load('my_est_params_pac_dss_ecuador');
end

guess_params_argentina = [0.0001   64.8784   -3.2567   25.7968];
guess_params_brazil = [0.0100  107.7532   -3.2013    0.8010];
guess_params_ecuador = [1e-5   21.82   -4.30   61.58];

%set target for moments from the data, and initial guesses for estimated parameters
if Argentina == 1
    target = [(4.77);(1.31);(3.81);1.78]; %Target moments from Argentina data
    guess_params = guess_params_argentina;
elseif Brazil == 1
    target = [4.64;1.1;1.65;0.1]; %Target moments from Brazil data
    guess_params = guess_params_brazil;
elseif Ecuador == 1
    target = [2.46;2.48;9.32;3.86]; %Target moments from Ecuador data
    guess_params = guess_params_ecuador; % optimal vector, distance = 2.36749
end

if exist('my_est_params','var')
    guess_params = my_est_params;
end

%Make monotonic change of the guess (ensures routine picks values within
%domain)
guess_params_opt = guess_params;
guess_params_opt(1) = log(guess_params_opt(1));
guess_params_opt(2) = (guess_params_opt(2))^0.5;


%Estimation function. It returns, for each targeted moment, the difference
%between target and simulation
myfun=@(x)my_smm_diff_function_dss(x,target,replications,innovations,model,params,M,approx1,T0,Ti,1);
options = optimoptions('lsqnonlin','Display','iter');

%First, compute simulated moments given the initial estimation guess, compare them
%with targets, and evaluate objective function (squared sum of differences)
sim_moments = target + myfun(guess_params_opt);
my_comp = [target sim_moments]
obj_fun = (sim_moments-target)'*(sim_moments-target)
% 
%Then, execute the whole estimation
%Lower and upper bounds for parameters
my_lb = [log(1e-5), 0, -Inf, -Inf];
my_ub = [log(0.01), Inf, Inf, Inf];

if perfom_estimate ==1
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
D1 = (my_est(4)) %DSS of debt

%Report full set of moments of interest
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1)
my_est_params = [phi_d phipar sigma_x D1]

%save the results
if save_estimate==1
    if Argentina == 1
        save_path = fullfile(my_root,'estimated parameters\my_est_params_pac_dss_argentina.mat');
    elseif Brazil == 1
        save_path = fullfile(my_root,'estimated parameters\my_est_params_pac_dss_brazil.mat');
    elseif Ecuador == 1
        save_path = fullfile(my_root,'estimated parameters\my_est_params_pac_dss_ecuador.mat');
    end
    save(save_path,'my_est_params');
end

if compute_errors==1
    fgru_pac_dss_errors;
end