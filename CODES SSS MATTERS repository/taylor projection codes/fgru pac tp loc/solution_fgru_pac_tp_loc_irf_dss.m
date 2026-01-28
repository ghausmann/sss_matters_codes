%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the DSS
% Solution method is Taylor Projection
% Calibration for Argentina, using DSS solution
% Computes consumption elasticity (Table 4 in the paper)
%---------------------------------------------------------------------------

clear;

%set the type of shock; you must not set more than one entry to 1 
irf_level = 0;
irf_prod = 0;
irf_vol = 1;
%set origin of estimation parameters (0 if full control, 1 if you load
%them)
origin_est = 1;
%Choose color for plots
my_color = 'b';
%Start fresh with the search path, and activate the required subfolders
my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
methods_path = genpath(fullfile(my_root,'taylor projection codes\solution_methods'));
model_path = genpath(fullfile(my_root,'taylor projection codes\model pac'));
addpath(methods_path);
addpath(model_path);
addpath(fullfile(my_root,'estimated parameters'));
addpath('files'); 
load('model')

%Calibration 
%Common parameters:
risk_factor = 0;
rho_eps_tb=0.95;
sigma_tb_bar=-8.06;     
rho_sigma_tb=0.94;
eta_tb=0.13;
nu=5;                   
etapar = 1001;
delta=0.014;            
alppha = 0.32;          
rho_x=0.95;             

%Country-specific parameters (Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);

if origin_est==1
    load('my_est_params_pac_dss_argentina'); %load estimated parameters
    phi_d = my_est_params(1); %controls Uzawa
    phipar = my_est_params(2); %controls capital adjustment costs
    sigma_x = log(exp(my_est_params(3))); %std. of innovations to productivity
    D_bar = my_est_params(4); %SSS debt
    
else
    phi_d = 0.2336;
    phipar = 28.5180;
    sigma_x = -2.8842;
    D_bar = 2.4150;
    
end

my_est_params
betta = 1/(1+(r_bar));   %discount factor

%DSS of auxiliary model
Hstar = ( (1-alppha)*(alppha/((r_bar)+(delta)))^(alppha/(1-alppha))  )^(1/(etapar-1));
Kstar = Hstar*(alppha/((r_bar)+(delta)))^(1/(1-alppha));
Istar = (delta)*Kstar;
Ystar = (Kstar^alppha)*Hstar^(1-alppha);
Cstar = D_bar/(1+(r_bar))-D_bar+Ystar-Istar;
C_bar = Cstar;
lambdastar = (Cstar)^(-nu);
phistar = lambdastar;
lKg = log(Kstar);

% Make a vector of parameter values:
params=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, etapar, phi_d, C_bar, risk_factor];

% Compute the cross moments:
n_e=5; % number of shocks.
M2f = eye(n_e);
M2f(1,1) = 0.00001; %This makes uncertainty shocks unexpected
M2f(2,2) = 0.00001;
[n_nodes,nodes,weights] = Monomials_1(n_e,M2f); % this quadrature function was written by Judd, Maliar, Maliar and Valero (2014).
nodes=nodes'; % transpose to n_e-by-n_nodes
M=get_moments(nodes,weights,model.order(2));
M.M2 = M2f;
%M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.

%----------------------------------------------------------------------
% Prepare an initial guess - in this case I use a perturbation solution
%----------------------------------------------------------------------

%DSS for states
nxss=[D_bar;log(Kstar);log(Istar);sigma_r_bar;sigma_tb_bar;0;0;0];
%DSS for controls
nyss=[log(Cstar);log(Hstar);(r_bar);log(phistar);log(Istar);D_bar];

% Compute the perturbation solution (keep the 4 outputs):
[derivs,stoch_pert,nonstoch_pert,model]=get_pert(model,params,M,nxss,nyss);

%-------------------------------------
% Solve the model by Taylor projection
%-------------------------------------

x0=nxss; % the initial approximation point
c0=nxss; % the center of the initial guess

% tolerance parameters for the Newton solver
tolX=1e-6;
tolF=1e-6;
maxiter=30;
%model.jacobian='exact'; % this is the default
model.jacobian='approximate'; % for large models try the approximate jacobian.
initial_guess=nonstoch_pert;
%initial_guess=stoch_pert;

% ---- Homotopy step -------
risk_factor_target = 1;
coeffs = initial_guess;
for h=0:1
    risk_factor1 = (1-h)*risk_factor + h*risk_factor_target;
    params(symparams=='risk_factor') = risk_factor1;
    [coeffs,model]=tpsolve(coeffs,x0,model,params,c0,nodes,weights,tolX,tolF,maxiter);
end

my_comp = [nonstoch_pert(1:20) stoch_pert(1:20) coeffs(1:20)];

%------------------------------------------------------------------
% Compute Taylor Projection solution
%------------------------------------------------------------------
    
[coeffs,model]=tpsolve(coeffs,x0,model,params,c0,nodes,weights,tolX,tolF,maxiter);
    
%---------------------------------
% Compute derivatives (for Table 4)
%---------------------------------

my_h = 1e-5;
x0 = nxss;
x1 = nxss; x1(4)=x1(4)+my_h;
g0_3=evalg(x0,coeffs,c0);
g1_3=evalg(x1,coeffs,c0);
%First number is the consumption elasticity
my_deriv_3 = (g1_3-g0_3)/my_h

