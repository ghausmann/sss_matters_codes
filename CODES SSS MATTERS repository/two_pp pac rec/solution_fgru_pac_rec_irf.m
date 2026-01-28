%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Auxiliary model (SSS recalibrated "rec" solution)
% Computes IRFs (for figures 3, 4, and 5)
%---------------------------------------------------------------------------

clear;

%set the type of shock; you must not set more than one entry to 1 
irf_level = 0;
irf_prod = 0;
irf_vol = 1;
%set origin of estimation parameters (0 if full control, 1 if you load
%them)
origin_est = 1;

approx0 = 2; %order of approximation for SSS
approx1 = 3; %order of approximation for solution
% Choose algorithm for the Sylvester equation
algo='vectorize'; % Simple vectorization.
%algo='dlyap'; % Hessenberg-Schur algorithm.
%algo='gensylv'; % Kamenik algorithm

%Choose color for plots
my_color = 'g:';

%Start fresh with the search path, and activate the required subfolders
my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
methods_path = genpath(fullfile(my_root,'sss_matters_methods'));
model_path = genpath(fullfile(my_root,'model pac sss'));
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
%nu=5;                  %Note that we estimate this one now
etapar= 1001;
delta=0.014;            
alppha = 0.32;          
rho_x=0.99;             %CALIBRATION CHANGE
%rho_x=0.95;            
rho_eps = 1;            
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
load('my_est_params_pac_rec'); %load estimated parameters
phi_d = my_est_params(1); %controls Uzawa
phipar = my_est_params(2); %controls capital adjustment costs
sigma_x = my_est_params(3); %std. of innovations to productivity
D1 = my_est_params(4); %SSS debt
nu = my_est_params(5);

my_est_params

D_bar = D1;
betta = 1/(1+(r_bar));   %discount factor

%Make vector of parameter values and matrices of innovations
routine_rec_tasks;
eps_ind = find(x=='eps'); %index of epsilon in state vector xt

%--------------------------------------------------------------------------
% Solving the model around the SSS 
%--------------------------------------------------------------------------

%Calibration of SSS (implied betta in model of interest, log(Ksss), excess
%depreciation in auxiliary model)
sss_sol = compute_betta_fgru_3p(D1,model,params,my_dss_params,M,eta,eps_ind,approx0,[-0.01 0.9*log(Kstar0) 0])
params(symparams=='psibetta') = sss_sol(1);
lK1 = sss_sol(2);
%Check that the SSS calibration delivers the targeted SSS
check_sss = compute_sss_fgru_3p(model,params,my_dss_params,M,eta,eps_ind,approx0,[D1 lK1 0])

%Compute aux. DSS and change parameter values as needed
[nxss,nyss,psi_k,C_bar] = my_dss_aux([D1 sss_sol(2:3)],my_dss_params)
params(symparams=='D_bar') = D1;
params(symparams=='psi_k') = psi_k;
params(symparams=='C_bar') = C_bar;
params(symparams=='psi_i') = psi_i;
% Compute the perturbation solution:
%DSS check
mycheck0 = double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
mycheck=sum(abs(mycheck0))

%Calculate derivatives of the policy functions abouts the DSS
derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx1,algo);

%------------------------------------------------------------------------------
% SSS results
%------------------------------------------------------------------------------

x0=nxss; % start at the steady state
derivs.hx(eps_ind,eps_ind) = 1;
x0(eps_ind) = 1; %impose the model of interest (epsilon=1)

%Evaluate decision rules for states and controls
mh_sss = dr_ht(derivs,nxss,approx0,(x0-nxss));
mg_sss = dr_gt(derivs,nyss,approx0,(x0-nxss));
y_sss = exp(mg_sss(y=='Y'));
c_sss = exp(mg_sss(y=='C'));
i_sss = exp(mg_sss(y=='If'));
D_sss = mg_sss(y=='Df');
r_sss = mg_sss(y=='r');

%SSS net exports ratio
nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss
%SSS net exports ratio (alternative way to compute it)
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) + 0.5*phi_d*(D_sss-D1)^2   )/y_sss
%Annualized external debt at SSS
nfa_sss = 100*D_sss/(12*y_sss)

%--------------------------------------------------------------------------
% Simulation of the model of interest (irfs)
%--------------------------------------------------------------------------
T = 100;
Ts = 9;

shocks = zeros(n_e,T); % a bunch of zeros
if irf_prod == 1
    shocks(5,Ts) = 1; %TFP shock
elseif irf_vol == 1
    shocks(1,Ts) = 1; %uncertainty shock
elseif irf_level == 1
    shocks(3,Ts) = 1; %interest rate level shock
end

%simulate the model
[yt,xt]=simul(x0,shocks,nyss,nxss,eta,derivs,approx1,0,model);

% in percent deviations from SSS (except possibly debt, just in deviations)
Ct = 100*(exp(yt(y=='C',:)) - c_sss)/c_sss;
Yt = 100*(exp(yt(y=='Y',:)) - y_sss)/y_sss;
Dt1 = (yt(y=='Df',:) - D_sss)/y_sss;
It1 = 100*(exp(yt(y=='If',:)) - i_sss)/i_sss;
rt = 100*(yt(y=='r',:) - r_bar);

%Compute Euler errors
betta_moi = betta*(((1+sss_sol(1))));
P2 = [betta_moi nu phi_d D_bar];
Sigma = M2f;
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_2(n_e,Sigma);
[~ ,mT] = size(xt);
lerrors =zeros(1,mT);
cum_errors = zeros(1,T-(Ts+1));
cum_time = Ts+1:T+1;

tic
for t=1:mT
    lerrors(t) =  log10(euler_errors_fgru_pac(P2,nxss,nyss,derivs,xt(:,t),eta,epsi_nodes,weight_nodes,approx1));
    
    if t>Ts
    cum_errors(t-Ts) = mean(lerrors(Ts+1:t));
    end
end
toc

errors_stats = ([mean(lerrors) max(lerrors)])

figure;
subplot(2,3,1);plot(Ct,my_color,'LineWidth',2);title('Consumption', 'FontSize', 12);
subplot(2,3,2);plot(Yt,my_color,'LineWidth',2);title('Output','FontSize',12);
subplot(2,3,3);plot(It1,my_color,'LineWidth',2);title('Investment','FontSize',12);
subplot(2,3,4);plot(Dt1,my_color,'LineWidth',2);title('Debt','FontSize',12);
subplot(2,3,5);plot(rt,my_color,'LineWidth',2);title('Interest rate','FontSize',12);
subplot(2,3,6);plot(cum_time,cum_errors,my_color,'LineWidth',2);title('Cumulative avrg. errors','FontSize',12);

