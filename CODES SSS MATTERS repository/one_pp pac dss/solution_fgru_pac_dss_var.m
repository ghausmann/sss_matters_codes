%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the DSS
% Computes simulated moments for variance decomposition (Table 2 in the
% paper)
%---------------------------------------------------------------------------

clear;
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

load('my_est_params_pac_dss_argentina'); %load estimated parameters
phi_d = my_est_params(1); %controls Uzawa
phipar = my_est_params(2); %controls capital adjustment costs
sigma_x = my_est_params(3); %std. of innovations to productivity
D1 = my_est_params(4); %SSS debt

my_est_params

D_bar = D1;
betta = 1/(1+(r_bar));   %discount factor

%Perform routine DSS and moment-setting tasks
routine_est_dss_tasks;

%--------------------------------------------------------------------------
% Estimation of the model
%--------------------------------------------------------------------------

%Generate matrix with shocks
T0 = 2000; %initial number of periods in simulations (before innovations occur) THIS IS INPUT FOR ESTIMATION AND MOMENT FUNCTIONS
Tt = 95;
Ti = 100;
T = Ti+Tt;
replications = 200;
sP1 = [rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb];
sP2 = [rho_eps_r sigma_r_bar rho_sigma_r eta_r r_bar 0];

guess_params = my_est_params;

%Make monotonic change of the guess (ensures routine picks values within
%domain)
guess_params_opt = guess_params;
guess_params_opt(1) = log(guess_params_opt(1));
guess_params_opt(2) = (guess_params_opt(2))^0.5;

my_est = guess_params_opt;

%all shocks
rng(1) %fix seed
[innovations,~] = my_innovations_var_decomp(T,replications,sP1,sP2,0,0,0);
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1);

var_y = yf(1);
var_c = yf(2)*var_y;
var_i = yf(3)*var_y;

my_var_decomp(:,1) = [var_y;var_c;var_i];

%TFP only
rng(1) %fix seed
[innovations,~] = my_innovations_var_decomp(T,replications,sP1,sP2,1,1,0);
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1);

var_y = yf(1);
var_c = yf(2)*var_y;
var_i = yf(3)*var_y;

my_var_decomp(:,2) = [var_y;var_c;var_i];

%w/o uncertainty
rng(1) %fix seed
[innovations,~] = my_innovations_var_decomp(T,replications,sP1,sP2,1,0,0);
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1);

var_y = yf(1);
var_c = yf(2)*var_y;
var_i = yf(3)*var_y;

my_var_decomp(:,3) = [var_y;var_c;var_i];

%rate level only
rng(1) %fix seed
[innovations,~] = my_innovations_var_decomp(T,replications,sP1,sP2,1,0,1);
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1);

var_y = yf(1);
var_c = yf(2)*var_y;
var_i = yf(3)*var_y;

my_var_decomp(:,4) = [var_y;var_c;var_i];

%w/o TFP
rng(1) %fix seed
[innovations,~] = my_innovations_var_decomp(T,replications,sP1,sP2,0,0,1);
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1);

var_y = yf(1);
var_c = yf(2)*var_y;
var_i = yf(3)*var_y;

my_var_decomp(:,5) = [var_y;var_c;var_i];

%uncertainty only
rng(1) %fix seed
[innovations,~] = my_innovations_var_decomp(T,replications,sP1,sP2,0,1,1);
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1);

var_y = yf(1);
var_c = yf(2)*var_y;
var_i = yf(3)*var_y;

my_var_decomp(:,6) = [var_y;var_c;var_i];

my_var_decomp