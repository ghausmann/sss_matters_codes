%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the DSS
% This script computes IRFs (for figures 2, 3, 5, 6, 7, 9 and 10 in paper)
%---------------------------------------------------------------------------

clear;

%set the type of shock; you must not set more than one entry to 1 
irf_level = 0;
irf_prod = 0;
irf_vol = 1;
%set the country; you must not set more than one entry to 1
Argentina = 1;
Brazil = 0;
Ecuador = 0;
%set origin of estimation parameters (0 if full control, 1 if you load
%them)
origin_est = 1;
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
    
    if origin_est==1
        load('my_est_params_pac_dss_argentina'); %load estimated parameters
        phi_d = my_est_params(1); %controls PAC
        phipar = my_est_params(2); %controls capital adjustment costs
        sigma_x = my_est_params(3); %std. of innovations to productivity
        D1 = my_est_params(4); %DSS debt
    else
        phi_d = 7e-5;
        phipar = 64.88;
        sigma_x = log(0.039);
        D1 = 25.8;
    end
elseif Brazil ==1
    rho_eps_r=0.95;
    sigma_r_bar=-6.97;
    rho_sigma_r=0.95;
    eta_r=0.28;
    r_bar=(0.007);
    
    if origin_est==1
        load('my_est_params_pac_dss_brazil'); %load estimated parameters
        phi_d = my_est_params(1); 
        phipar = my_est_params(2); 
        sigma_x = my_est_params(3); 
        D1 = my_est_params(4); 
    else
        phi_d = 0.01;
        phipar = 107.75;
        sigma_x = log(0.04);
        D1 = 0.8;
    end
elseif Ecuador == 1
    rho_eps_r=0.95;
    sigma_r_bar=-6.06;
    rho_sigma_r=0.96; 
    eta_r=0.35;
    r_bar=0.011;
        
    if origin_est==1
        load('my_est_params_pac_dss_ecuador'); %load estimated parameters
        phi_d = my_est_params(1); 
        phipar = my_est_params(2); 
        sigma_x = my_est_params(3); 
        D1 = my_est_params(4); 
    else
        phi_d = 1e-5;
        phipar = 19.34;
        sigma_x = log(0.014);
        D1 = 58.63;
    end

end

D_bar = D1;
betta = 1/(1+(r_bar));   

%Compute the perturbation solution
routine_sol_dss_tasks;

%------------------------------------------------------------------------------
% simulate the model for Tsss periods (FIRST RUN to converge to an implied SSS)
%------------------------------------------------------------------------------
x0=nxss; % start at the DSS
Tsss = 2000;
shocks0 = zeros(n_e,Tsss); % a bunch of zeros
[myt,mxt]=simul(x0,shocks0,nyss,nxss,eta,derivs,approx1,0,model);
xsss = mxt(:,Tsss+1:end);
ysss = myt(:,Tsss+1:end);

%Implied SSS values
y_sss = exp(ysss(5));
c_sss = exp(ysss(1));
i_sss = exp(ysss(7));
D_sss = ysss(8);
r_sss = ysss(3);

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
%Deterministic simulation of annualized external debt
nfa = 100*myt(8,:)./(12*exp(myt(5,:)));
figure;plot(nfa(1:600));xlabel('Time t');

%--------------------------------------------------------------------------
% Simulation to compute irfs
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

[yt,xt]=simul(xsss,shocks,nyss,nxss,eta,derivs,approx1,0,model);

% in percent deviations from SSS (except debt, in deviations-to-SSS_output ratio)
Ct = 100*(exp(yt(y=='C',:)) - c_sss)/c_sss;
Yt = 100*(exp(yt(y=='Y',:)) - y_sss)/y_sss;
Dt1 = (yt(y=='Df',:) - D_sss)/y_sss;
It1 = 100*(exp(yt(y=='If',:)) - i_sss)/i_sss;
rt = 100*(yt(y=='r',:) - r_bar);

%Compute Cumulative Euler errors
P2 = [betta nu phi_d D_bar];
Sigma = M2f;
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_1(n_e,Sigma);
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

figure;
subplot(2,3,1);plot(Ct,my_color,'LineWidth',2);title('Consumption', 'FontSize', 12);
subplot(2,3,2);plot(Yt,my_color,'LineWidth',2);title('Output','FontSize',12);
subplot(2,3,3);plot(It1,my_color,'LineWidth',2);title('Investment','FontSize',12);
subplot(2,3,4);plot(Dt1,my_color,'LineWidth',2);title('Debt','FontSize',12);
subplot(2,3,5);plot(rt,my_color,'LineWidth',2);title('Interest rate','FontSize',12);
subplot(2,3,6);plot(cum_time,cum_errors,my_color,'LineWidth',2);title('Cumulative avrg. errors','FontSize',12);

