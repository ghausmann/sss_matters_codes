%---------------------------------------------------------------------------
% FGRU Small open-economy model
% UZAWA model, approximation around the DSS
% Estimation of parameters (Table 5 in the paper)
%---------------------------------------------------------------------------

clear;

%Choose whether to perform a full estimation (1), or just report moments (0)
perfom_estimate = 0;
%Choose whether to save estimation results
save_estimate = 0;
%Choose whether to also compute Euler errors and deterministic simulation
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
betta = 1/(1+(r_bar));   %discount factor

%To-be estimated parameters.
load('my_est_params_uza_dss'); %load estimated parameters
phi_d = my_est_params(1); %controls Uzawa
phipar = my_est_params(2); %controls capital adjustment costs
sigma_x = my_est_params(3); %std. of innovations to productivity
D_bar = my_est_params(4); %SSS debt
D1 = D_bar;

%DSS of the model
Hstar = ( (1-alppha)*(alppha/((r_bar)+(delta)))^(alppha/(1-alppha))  )^(1/(etapar-1));
Kstar = Hstar*(alppha/((r_bar)+(delta)))^(1/(1-alppha));
Istar = delta*Kstar;
Ystar = (Kstar^alppha)*Hstar^(1-alppha);
Cstar = D1/(1+(r_bar))-D1+Ystar-Istar;
C_bar = Cstar;
lambdastar = (Cstar)^(-nu);
phistar = lambdastar;

nxss=[D1;log(Kstar);log(Istar);sigma_r_bar;sigma_tb_bar;0;0;0;0;0;0];
nyss=[log(Cstar);log(Hstar);(r_bar);log(lambdastar);log(Ystar);log(phistar);log(Istar);D1];


% Make a vector of parameter values:
params=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, rho_eps, etapar, phi_d, C_bar, psibetta, psi_k, psi_i];

% Make eta
% Define eta as in Schmitt-Grohe and Uribe (2004):
eta=[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;eta_r 0 0 0 0;0 eta_tb 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 exp(sigma_x);0 0 0 0 0];

% Compute the cross moments:
n_e=size(eta,2); % number of shocks.
M2f = eye(n_e);
M2f(1,1) = 0.00001; %This makes uncertainty shocks unexpected (only relevant for fifth-order solutions)
M2f(2,2) = 0.00001;
[n_nodes,nodes,weights] = Monomials_2(n_e,M2f); % this quadrature function was written by Judd, Maliar, Maliar and Valero (2014).
nodes=nodes'; % transpose to n_e-by-n_nodes
M=get_moments(nodes,weights,approx1);

%--------------------------------------------------------------------------
% Estimation of the model
%--------------------------------------------------------------------------

%Generate matrix with shocks
T0 = 2000; %initial number of periods in simulations (before innovations occur) THIS IS INPUT FOR ESTIMATION AND MOMENT FUNCTIONS
Tt = 95; %number of periods in each stochastic simulation
Ti = 100; %burn-in number of periods
T = Ti+Tt;
replications = 200;
sP1 = [rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb];
sP2 = [rho_eps_r sigma_r_bar rho_sigma_r eta_r r_bar 0];
[innovations,~] = my_innovations(T,replications,sP1,sP2);

%set target for moments from the data, and initial guesses for estimated parameters
target = [(4.77);(1.31);(3.81);1.78]; %Target moments from Argentina data
guess_params = my_est_params

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
    %Uncomment this instead to use as estimated parameters the initial guess
    my_est = guess_params_opt;
end

%estimated parameters
phi_d = exp(my_est(1)); %controls PAC
phipar=my_est(2)^2; %capital adjustment costs
sigma_x = (my_est(3)); %std. innovations to productivity
D1 = (my_est(4)) %DSS of debt

%Report full set of moments of interest
yf = my_smm_moments_dss_compare(my_est,replications,innovations,model,params,M,approx1,T0,Ti,1)

%save the results
my_est_params = [phi_d phipar sigma_x D1]
if save_estimate==1
    save_path = fullfile(my_root,'estimated parameters\my_est_params_uza_dss.mat');
    save(save_path,'my_est_params');
end

if compute_errors==1
    fgru_uza_dss_errors;
end