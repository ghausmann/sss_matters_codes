%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Auxiliary model (SSS fixed), approximation around the SSS
% Computes IRFs (for figures 2, 6, 7 and 11 in paper)
%---------------------------------------------------------------------------

clear;

%set the type of shock; you must not set more than one entry to 1 
irf_level = 0;
irf_prod = 0;
irf_vol = 1;
%set origin of estimation parameters (0 if full control, 1 if you load
%them)
origin_est = 1;
%Choose if you want to activate a robustnetss check (default is not) 
check_phi_d = 0; %fixed PAC parameter
check_nu = 0; %log utility

%For a fifth-order solution, set approx0 = 4 and approx1 = 5 for speed
approx0 = 2; %order of approximation for SSS
approx1 = 3; %order of approximation for solution
% Choose algorithm for the Sylvester equation
algo='vectorize'; % Simple vectorization.
%algo='dlyap'; % Hessenberg-Schur algorithm.
%algo='gensylv'; % Kamenik algorithm

%Choose color for plots
my_color = 'r--';

%Start fresh with the search path, and activate the required subfolders
my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
methods_path = genpath(fullfile(my_root,'sss_matters_methods'));
model_path = genpath(fullfile(my_root,'model pac loc'));
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
%Auxiliary parameters
psi_k = 0;
psi_i= 0; 
Da = 0;
%Set auxiliary PAC parameter
if check_nu==1
    psi_d = 0.0014;
elseif check_phi_d==1
    psi_d=6e-4;
else
    psi_d=4e-4;
end

%Country-specific parameters (Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);
if origin_est==1
    if check_nu==1
        load('my_est_params_pac_dss_argentina_nu');
    elseif check_phi_d==1
        load('my_est_params_pac_dss_argentina_phi_d');
    else
        load('my_est_params_pac_dss_argentina'); %load estimated parameters
    end
    
    phi_d = my_est_params(1); %controls Uzawa
    phipar = my_est_params(2); %controls capital adjustment costs
    sigma_x = my_est_params(3); %std. of innovations to productivity
    D1 = my_est_params(4); %SSS debt
else
    phi_d = 0.0002;
    phipar = 64.5975;
    sigma_x = -2.8768;
    D1 = 21.6190;
end

my_est_params

D_bar = D1;
betta = 1/(1+(r_bar));   %discount factor

%Make vector of parameter values and matrices of innovations
routine_sss_loc_tasks;
eps_ind = find(x=='eps'); %index of epsilon in state vector xt

%--------------------------------------------------------------------------
% Compute SSS and dynamic solution around it
%--------------------------------------------------------------------------

%Solve for the SSS
t =clock;
sss_sol = compute_sss_fgru_3_locp_num(model,params,my_dss_params,M,eta,eps_ind,approx0,[1 log(Kstar0) 0])
my_time1 = etime(clock,t)

%Recompute DSS of auxiliary model
[nxss,nyss,psi_k] = my_num_dss(sss_sol,my_dss_params);

params(symparams=='Da') = sss_sol(1);
params(symparams=='psi_k') = psi_k;
params(symparams=='psi_i') = sss_sol(3);

% Compute the perturbation solution:
%DSS check
mycheck0 = double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
mycheck=sum(abs(mycheck0))
%Calculate derivatives of the policy functions abouts the DSS
derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx1,algo);

%--------------------------------------------------------------------------
% SSS results
%--------------------------------------------------------------------------

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
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) + 0.5*phi_d*(D_sss-D_bar)^2   )/y_sss
%Punishment for not holding a D1 exogenous debt level
crazy_punish = 100*(  0.5*phi_d*(D_sss-D_bar)^2   )/y_sss
%Annualized external debt at DSS
nfa_dss = 100*nyss(8)/(12*exp(nyss(5)))
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
betta_moi = (1/(1+r_bar));
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

%Metric to assess accuracy of auxiliary model under the choice of psi_d
%is expected portfolio euler errors conditional on being at the SSS
[NN ,~] = size(epsi_nodes);
my_errors =zeros(1,NN);
my_x0 = x0 + eta*epsi_nodes';

for nn=1:NN
    
    my_errors(nn) = euler_errors_fgru_pac(P2,nxss,nyss,derivs,my_x0(:,nn),eta,epsi_nodes,weight_nodes,approx1);
    
end

sss_error = lerrors(1)
my_weighted_errors = weight_nodes'*my_errors'
l_weighted_errors = log10(my_weighted_errors)


figure;
subplot(2,3,1);plot(Ct,my_color,'LineWidth',2);title('Consumption', 'FontSize', 12);
subplot(2,3,2);plot(Yt,my_color,'LineWidth',2);title('Output','FontSize',12);
subplot(2,3,3);plot(It1,my_color,'LineWidth',2);title('Investment','FontSize',12);
subplot(2,3,4);plot(Dt1,my_color,'LineWidth',2);title('Debt','FontSize',12);
subplot(2,3,5);plot(rt,my_color,'LineWidth',2);title('Interest rate','FontSize',12);
subplot(2,3,6);plot(cum_time,cum_errors,my_color,'LineWidth',2);title('Cumulative avrg. errors','FontSize',12);

