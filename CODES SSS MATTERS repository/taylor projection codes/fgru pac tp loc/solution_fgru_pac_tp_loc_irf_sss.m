%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the SSS
% Solution method is Taylor Projection
% Calibration for Argentina, using DSS solution
% Computes IRFs (figures 6 and 7 in paper) and elasticities (Table 4) 
% NOTE: This version of the codes slightly improves the Euler errors reported in the
% paper.
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
% Compute SSS by iteration
%------------------------------------------------------------------
T = 2000;
shocks = zeros(5,T+1); % a bunch of zeros

nguess = [D_bar;log(Kstar);log(Istar)]
x0 = nxss;
x0(1:3) = nguess(1:3);

my_dist=1;

while my_dist>=1e-5
    
    [coeffs,model]=tpsolve(coeffs,x0,model,params,c0,nodes,weights,tolX,tolF,maxiter);
    
    % preallocate
    x_simul0=zeros(model.n_x,T+1);
    y_simul0=zeros(model.n_y,T);
    x_simul0(:,1)=x0;
    
   
    for t=1:T
        xt=x_simul0(:,t);
        up0=shocks(:,t+1);
               
        yt=evalg(xt,coeffs,c0);
        
        y_simul0(:,t)=yt;
        x_simul0(:,t+1)=evalPhi(xt,yt,up0,params);
             
    end
    
    my_dist = norm(x_simul0(1:3,end)-x0(1:3))
    x0 = x_simul0(:,end);
        
end

SSS_Taylor = x0

%------------------------------------------------------------------
% Compute the residual function and the model variables at point x0
%------------------------------------------------------------------

[R_fun0,g_fun0,Phi_fun0,auxvars0]=residual(coeffs,x0,params,c0,nodes,weights);
% compute the function g(x) at x0
y0=evalg(x0,coeffs,c0)
% 
% compute the function Phi(x,y,epsp) at x0, y0 and epsp0
up0=[0;0;0;0;0];
xp0=evalPhi(x0,y0,up0,params)
diff_x0 = xp0-x0

%---------------------------------
% simulate the model for T periods
%---------------------------------
T = 101;
Ts = 10;

shocks = zeros(5,T+1); % a bunch of zeros
if irf_prod == 1
    shocks(5,Ts) = 1; %TFP shock
elseif irf_vol == 1
    shocks(1,Ts) = 1; %uncertainty shock
elseif irf_level == 1
    shocks(3,Ts) = 1; %interest rate level shock
end

% preallocate
x_simul=zeros(model.n_x,T+1);
y_simul=zeros(model.n_y,T);
R_simul=zeros(model.n_y,T);
lYt = zeros(1,T);
lerrors = zeros(1,T);
cum_errors = zeros(1,T-(Ts+1));
cum_time = Ts:T;

x_simul(:,1)=x0;

%option=1; % compute only simulated variables
option=2; % compute model residuals

for t=1:T
    xt=x_simul(:,t);
    up0=shocks(:,t+1);

    % Option 1 - compute only the simulated variables
    if option==1
        yt=evalg(xt,coeffs,c0);
        
        y_simul(:,t)=yt;
        x_simul(:,t+1)=evalPhi(xt,yt,up0,params);
    else
    % Option 2 - compute also model residuals
        [Rt,yt,~,auxt]=residual(coeffs,xt,params,c0,nodes,weights);
        
        y_simul(:,t)=yt;
        x_simul(:,t+1)=evalPhi(xt,yt,up0,params);
        R_simul(:,t)=Rt;
        lerrors(t)=log10(abs(R_simul(1,t)));
        if t>Ts-1
            cum_errors(t-(Ts-1)) = mean(lerrors((Ts-1)+1:t));
        end
        
        lYt(t) = auxt(7,1);
    end
end

% in percent deviations from SSS (except possibly debt, just in deviations)
Ct = 100*(exp(y_simul(1,:)) - exp(y_simul(1,1)))/exp(y_simul(1,1));
Yt = 100*(lYt - lYt(1))/lYt(1);
Dt1 = (y_simul(6,:) - y_simul(6,1))/lYt(1);
It1 = 100*(exp(y_simul(5,:)) - exp(y_simul(5,1)))/exp(y_simul(5,1));
rt = 100*(y_simul(3,:) - r_bar);

y_sss = lYt(1);
c_sss = exp(y_simul(1,2));
i_sss = exp(y_simul(5,2));

D_sss = y_simul(6,2);
r_sss = y_simul(3,2);

nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) + 0.5*phi_d*(D_sss-D_bar)^2   )/y_sss


figure;
subplot(2,3,1);plot(Ct,my_color,'LineWidth',2);title('Consumption');
subplot(2,3,2);plot(Yt,my_color,'LineWidth',2);title('Output');
subplot(2,3,3);plot(It1,my_color,'LineWidth',2);title('Investment');
subplot(2,3,4);plot(Dt1,my_color,'LineWidth',2);title('Debt');
subplot(2,3,5);plot(rt,my_color,'LineWidth',2);title('Interest rate');
subplot(2,3,6);plot(cum_time,cum_errors,my_color,'LineWidth',2);title('Cumulative avrg. errors','FontSize',12);

%---------------------------------
% Compute derivatives (for Table 4)
%---------------------------------
my_h = 1e-5;
x1 = x0; x1(4)=x1(4)+my_h;
g0_3=evalg(x0,coeffs,c0);
g1_3=evalg(x1,coeffs,c0);
%First number is the consumption elasticity
my_deriv_taylor_sss = (g1_3-g0_3)/my_h



