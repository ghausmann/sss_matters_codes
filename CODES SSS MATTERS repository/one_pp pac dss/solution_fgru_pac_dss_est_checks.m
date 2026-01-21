%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the DSS
% Estimation of parameters for robustness checks (
%---------------------------------------------------------------------------

clear;

%set the type of check 
check_phi_d = 1; %fixed PAC parameter
check_nu = 0; %log utility

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

%Country-specific parameters (for Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);

betta = 1/(1+(r_bar));   %discount factor

%To-be estimated parameters. Here just initialize them with some numbers
phi_d = 4.2e-4; %controls PAC
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

target = [(4.77);(1.31);(3.81);1.78]; %Target moments from Argentina data

%load previously estimated parameters (if available)
if check_nu==1
    load('my_est_params_pac_dss_argentina_nu');
elseif check_phi_d==1
    load('my_est_params_pac_dss_argentina_phi_d');
else
    load('my_est_params_pac_dss_argentina');
end

if exist('my_est_params','var')
    guess_params = my_est_params;
else
    guess_params = [0.0001   64.8784   -3.2567   25.7968];
end

if check_phi_d==1
   guess_params(1) = phi_d;
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
if check_phi_d==1
    my_lb = [log(phi_d), 0, -Inf, -Inf];%Fix phi_d in routine
    my_ub = [log(phi_d), Inf, Inf, Inf];
else
    my_lb = [log(1e-5), 0, -Inf, -Inf];
    my_ub = [log(0.01), Inf, Inf, Inf];
end

if perfom_estimate ==1
    t =clock;
    my_est = lsqnonlin(myfun,guess_params_opt,my_lb,my_ub,options)
    my_time1 = etime(clock,t)
else
    %Uncomment this instead to use as estimated parameters the initial guess
    my_est = guess_params_opt;
end

%estimated parameters
phi_d = exp(my_est(1)); %controls PAC
phipar=my_est(2)^2; %capital adjustment costs
sigma_x = (my_est(3)); %std. innovations to productivity
D1 = (my_est(4)) %DSS of debt

%Report full set of moments of interest (would be nice to have a table
%here)
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1)
%ym = my_smm_moments_dss_compare_mod(my_est,replications,innovations,model,params,M)
%save the results
my_est_params = [phi_d phipar sigma_x D1]

if save_estimate==1
    if check_nu==1
        save_path = fullfile(my_root,'estimated parameters\my_est_params_pac_dss_argentina_nu.mat');
    elseif check_phi_d==1
        save_path = fullfile(my_root,'estimated parameters\my_est_params_pac_dss_argentina_phi_d.mat');
    end
    save(save_path,'my_est_params');
end

if compute_errors==1
    fgru_pac_dss_errors;
end