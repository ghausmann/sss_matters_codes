%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the SSS
% Computes IRFs (for figures 2, 3, and 5 in paper, excluding Euler errors)
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
methods_path = genpath(fullfile(my_root,'sss_matters_methods\common functions'));
funs_path = genpath(fullfile(my_root,'dynare codes\model pac loc\functions'));
addpath(methods_path);
addpath(funs_path);
addpath(fullfile(my_root,'estimated parameters'));
addpath('C:\dynare\5.5\matlab'); %change path accordingly
% Load the model:
load('my_pac_loc_model.mat');

if origin_est==1
    %load estimated parameters
    load('my_est_params_pac_dss_argentina');     
    phi_d = my_est_params(1); %controls Uzawa
    phipar = my_est_params(2); %controls capital adjustment costs
    sigma_x = my_est_params(3); %std. of innovations to productivity
    D_bar = my_est_params(4); %SSS debt
else
    phi_d = 0.0001;
    phipar = 65;
    sigma_x = -3.2567;
    D_bar = 26;
    
end

%DSS of ORIGINAL model (as implied by DSS calibration)
Hstar0 = ( (1-alppha)*(alppha/((r_bar)+(delta)))^(alppha/(1-alppha))  )^(1/(etapar-1));
Kstar0 = Hstar0*(alppha/((r_bar)+(delta)))^(1/(1-alppha));
Istar0 = (delta)*Kstar0;
Ystar0 = (Kstar0^alppha)*Hstar0^(1-alppha);
Cstar0 = D_bar/(1+(r_bar))-D_bar+Ystar0-Istar0
C_bar0 = Cstar0;
my_dss_params = [alppha r_bar delta etapar nu D_bar phi_d C_bar0 sigma_r_bar sigma_tb_bar];
M_.params(strcmp(M_.param_names,'C_bar')==1)    = C_bar0;
M_.params(strcmp(M_.param_names,'phi_d')==1)    = phi_d;
M_.params(strcmp(M_.param_names,'phipar')==1)   = phipar;
M_.params(strcmp(M_.param_names,'sigma_x')==1)  = sigma_x;
M_.params(strcmp(M_.param_names,'D_bar')==1)    = D_bar;

%Collect indices of variables in vectors (useful to implement the SSS
%algorithm)
eps_ind = 9; %index of perturbation variable WITHIN the states
st_ini_ind = find(strcmp(oo_.var_list,'sigma_r')); %index of first state variable
st_end_ind = find(strcmp(oo_.var_list,'eps')); %index of last state variable
D_ind = find(strcmp(oo_.var_list,'D')); %index of debt state
K_ind = find(strcmp(oo_.var_list,'K')); %index of debt state
I_ind = find(strcmp(oo_.var_list,'I')); %index of debt state
my_ind = [eps_ind;st_ini_ind;st_end_ind;D_ind;K_ind;I_ind];
%Number of exogenous shocks
nshocks = M_.exo_nbr; 

%--------------------------------------------------------------------------
% Compute SSS and dynamic solution around it
%--------------------------------------------------------------------------

%Recompute DSS
[yss,~] = my_num_dss([D_bar log(Kstar0) 0],my_dss_params);
oo_.steady_state = yss;

%Calculate derivatives of the policy functions abouts the DSS
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%--------------------------------------------------------------------------
% Simulation of model of interest (irfs)
%--------------------------------------------------------------------------
T0 = 2000; %larger number of periods to converge to an initial SSS
T = 101;
Ts = 10-1;

% start at the DSS
ysss0 = mdr.ys;

innovations = zeros(5,(T0 + (T-1)));

if irf_prod == 1
    innovations(5,T0+Ts) = 1; %TFP shock
elseif irf_vol == 1
    innovations(1,T0+Ts) = 1; %uncertainty shock
elseif irf_level == 1
    innovations(3,T0+Ts) = 1; %interest rate level shock
    
end

myt =simult_(M_,options_,ysss0,mdr,innovations',options_.order);
myt = myt(:,T0+1:end);
ysss = myt(:,1); %implied SSS after the long deterministic simulation

% % in percent deviations from SSS
idx0 = find(strcmp(oo_.var_list,'C'));lCt = myt(oo_.dr.order_var == idx0,:);c_sss = exp(ysss(oo_.dr.order_var == idx0));
idx0 = find(strcmp(oo_.var_list,'Y'));lYt = myt(oo_.dr.order_var == idx0,:);y_sss = exp(ysss(oo_.dr.order_var == idx0));
idx0 = find(strcmp(oo_.var_list,'Dp'));Dpt = myt(oo_.dr.order_var == idx0,:);D_sss = (ysss(oo_.dr.order_var == idx0));
idx0 = find(strcmp(oo_.var_list,'Ip'));lIpt = myt(oo_.dr.order_var == idx0,:);i_sss = exp(ysss(oo_.dr.order_var == idx0));
idx0 = find(strcmp(oo_.var_list,'r'));rt_0 = myt(oo_.dr.order_var == idx0,:);

% in percent deviations from SSS (except debt, just in deviations over sss output)
Ct = 100*(exp(lCt) - c_sss)/c_sss;
Yt = 100*(exp(lYt) - y_sss)/y_sss;
Dt1 = (Dpt - D_sss)/y_sss;
It1 = 100*(exp(lIpt) - i_sss)/i_sss;
rt = 100*(rt_0 - r_bar);


figure;
subplot(2,3,1);plot(Ct,my_color,'LineWidth',2);title('Consumption', 'FontSize', 12);
subplot(2,3,2);plot(Yt,my_color,'LineWidth',2);title('Output','FontSize',12);
subplot(2,3,3);plot(It1,my_color,'LineWidth',2);title('Investment','FontSize',12);
subplot(2,3,4);plot(Dt1,my_color,'LineWidth',2);title('Debt','FontSize',12);
subplot(2,3,5);plot(rt,my_color,'LineWidth',2);title('Interest rate','FontSize',12);

