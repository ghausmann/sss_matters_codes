%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC, approximation around the DSS
% This script replicates Figure 1, panels (a) and (b) in the paper
%---------------------------------------------------------------------------

clear;
%Choose order of approximation
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

%Compute the perturbation solution
routine_sol_dss_tasks;

%------------------------------------------------------------------------------
% simulate the model for Tsss periods (FIRST RUN to converge to an implied SSS)
%------------------------------------------------------------------------------
x0=nxss; % start at the steady state
Tsss = 2000;
shocks0 = zeros(n_e,Tsss); % a bunch of zeros
[myt,mxt]=simul(x0,shocks0,nyss,nxss,eta,derivs,approx1,0,model);
xsss = mxt(:,end);
ysss = myt(:,end);

if norm(mxt(:,end)-mxt(:,end-1))>1e-8
    disp('WARNING! Convergence to SSS incomplete')
end

%Implied SSS values
y_sss = exp(ysss(5));
c_sss = exp(ysss(1));
i_sss = exp(ysss(7));
D_sss = ysss(8);
r_sss = ysss(3);

%SSS net exports ratio
nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss
%SSS net exports ratio (alternative way to compute it)
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) + 0.5*phi_d*(D_sss-D1)^2   )/y_sss
%Punishment for not holding a D1 exogenous debt level
crazy_punish = 100*(  0.5*phi_d*(D_sss-D1)^2   )/y_sss
%Annualized external debt at DSS
nfa_dss = 100*nyss(8)/(12*exp(nyss(5)))
%Annualized external debt at SSS
nfa_sss = 100*D_sss/(12*y_sss)
%Deterministic simulation of annualized external debt
nfa = 100*myt(8,:)./(12*exp(myt(5,:)));
figure;plot(nfa(1:600));xlabel('Time t');

%---------------------------------------------------------------------------
%Loop
%---------------------------------------------------------------------------

n_a = 25; %number of elements in the grids
Dmin = D_sss;
Dmax = D_bar;
D_grid = [Dmin Dmin + ((1:(n_a-1))/(n_a-1)) * (D_bar - Dmin)];
v_nx_y_dss = zeros(1,25);
v_nx_y_sss = zeros(1,25);
v_nfa_dss = zeros(1,25);

for n = 1:n_a
    
    D_n = D_grid(n);
    Cstar_n = D_n/(1+(r_bar))-D_n+Ystar-Istar;
    C_bar_n = Cstar_n;
    lambdastar_n = (Cstar_n)^(-nu);
    
    params(symparams=="D_bar")=D_n;
    params(symparams=="C_bar")=C_bar_n;
    nxss_n=[D_n;log(Kstar);log(Istar);sigma_r_bar;sigma_tb_bar;0;0;0;0;0;0];
    nyss_n=[log(Cstar_n);log(Hstar);(r_bar);log(lambdastar_n);log(Ystar);log(phistar);log(Istar);D_n];

    %Calculate derivatives of the policy functions abouts the DSS
    derivs_n=solve_dsge(model,params,M,eta,nxss_n,nyss_n,approx1,algo);

    x0_n=nxss_n; % start at the steady state
    [myt_n,mxt_n]=simul(x0_n,shocks0,nyss_n,nxss_n,eta,derivs_n,approx1,0,model);
    if norm(mxt_n(:,end)-mxt_n(:,end-1))>1e-5
        fprintf('WARNING! Convergence to SSS incomplete in %d.\n', n);
    end
    ysss_n = myt_n(:,Tsss+1:end);

    y_sss_n = exp(ysss_n(5));
    c_sss_n = exp(ysss_n(1));
    i_sss_n = exp(ysss_n(7));
    
    
    v_nx_y_dss(n) = 100*(Ystar - Cstar_n - Istar)/Ystar;
    v_nx_y_sss(n) = 100*(y_sss_n - c_sss_n - i_sss_n)/y_sss_n;
    v_nfa_dss(n) = 100*nyss_n(8)/(12*exp(nyss_n(5)));

end

figure;plot(v_nfa_dss,v_nx_y_sss);xlabel('nfa');title('dss nfa vs sss net exports')


