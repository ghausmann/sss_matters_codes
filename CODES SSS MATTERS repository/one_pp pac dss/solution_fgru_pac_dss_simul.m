%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the SSS
% This script is just to simulate the exogenous states, so we can feed
% Taylor projection's simulations with them. 
%---------------------------------------------------------------------------

clear;

%Choose order of approximation for solution
approx1 = 3; 
% Choose algorithm for the Sylvester equation
algo='vectorize'; % Simple vectorization.
%algo='dlyap'; % Hessenberg-Schur algorithm.
%algo='gensylv'; % Kamenik algorithm

%Choose color for plots
my_color = 'b';

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

%Argentina-specific parameters:
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);

%load previously estimated parameters
load('my_est_params_pac_dss_argentina');

phi_d = my_est_params(1); %controls PAC
phipar = my_est_params(2); %controls capital adjustment costs
sigma_x = my_est_params(3); %std. of innovations to productivity
D1 = my_est_params(4); %SSS debt
   
my_est_params

D_bar = D1;
betta = 1/(1+(r_bar));   %discount factor

%Compute the perturbation solution
routine_sol_tasks;

%------------------------------------------------------------------------------
% simulate the model for Tsss periods (FIRST RUN to converge to an implied SSS)
%------------------------------------------------------------------------------
x0=nxss; % start at the steady state
Tsss = 2000;
shocks0 = zeros(n_e,Tsss); % a bunch of zeros
[myt,mxt]=simul(x0,shocks0,nyss,nxss,eta,derivs,approx1,0,model);
xsss = mxt(:,Tsss+1:end);
ysss = myt(:,Tsss+1:end);

y_sss = exp(ysss(5));
c_sss = exp(ysss(1));
i_sss = exp(ysss(7));

D_sss = ysss(8);
r_sss = ysss(3);
nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss

%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------

Tt = 95;
Ti = 100;
T = Ti+Tt;
replications =200;
sP1 = [rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb];
sP2 = [rho_eps_r sigma_r_bar rho_sigma_r eta_r r_bar 0];
[innovations,~] = my_innovations(T,replications,sP1,sP2);
eps_r_t_1 = zeros(replications,T+1);
eps_tb_t_1 = zeros(replications,T+1);


for t=1:replications
    
    shocks = innovations(:,:,t);
    %shocks = zeros(5,T);
    [yt,xt]=simul(xsss,shocks,nyss,nxss,eta,derivs,approx1,0,model);
        eps_r_t_1(t,:) = xt(6,:);
    eps_tb_t_1(t,:) = xt(7,:);
    
end

save('my_sim.mat','eps_r_t_1','eps_tb_t_1');