%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC, approximation around the DSS
% Computes consumption elasticities about the DSS (Table 4 in the paper)
%---------------------------------------------------------------------------

clear;
approx1 = 5; %Order of approximation. Needs to be 5th here
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

%Compute the perturbation solution
routine_sol_dss_tasks;

%------------------------------------------------------------------------------
% Compute derivatives
%------------------------------------------------------------------------------
my_h = 1e-5;
x0 = nxss;
x1 = nxss; x1(4)=x1(4)+my_h;

g0_3 = dr_gt(derivs,nyss,3,(x0-nxss));
g1_3 = dr_gt(derivs,nyss,3,(x1-nxss));
my_deriv_3 = (g1_3-g0_3)/my_h;

g0_4 = dr_gt(derivs,nyss,4,(x0-nxss));
g1_4 = dr_gt(derivs,nyss,4,(x1-nxss));
my_deriv_4 = (g1_4-g0_4)/my_h;

g0_5 = dr_gt(derivs,nyss,5,(x0-nxss));
g1_5 = dr_gt(derivs,nyss,5,(x1-nxss));
my_deriv_5 = (g1_5-g0_5)/my_h;

%First row is consumption elasticities
compare_derivs_standard = [my_deriv_3 my_deriv_4 my_deriv_5]
save('my_derivs_pac_dss','compare_derivs_standard')
