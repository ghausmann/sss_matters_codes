%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, approximation around the SSS
% Solution method is Taylor Projection
% Calibration for Argentina, using DSS solution
% Computes simulated moments (Table 3 in the paper)
% NOTE: This version of the codes slightly improves the Euler errors reported in the
% paper.
%---------------------------------------------------------------------------

clear;rng(1)

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
load('my_sim.mat');

%set origin of estimation parameters (0 if full control, 1 if you load
%them)
origin_est = 1;

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

%------------------------------------------------------------------
% Compute the residual function and the model variables at point x0
%------------------------------------------------------------------

[R_fun0,g_fun0,Phi_fun0,auxvars0]=residual(coeffs,x0,params,c0,nodes,weights);
% compute the function g(x) at x0
y0=evalg(x0,coeffs,c0)

y_sss = auxvars0(7);
c_sss = exp(y0(1));
i_sss = exp(y0(5));
nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss;

% 
% compute the function Phi(x,y,epsp) at x0, y0 and epsp0
up0=[0;0;0;0;0];
xp0=evalPhi(x0,y0,up0,params)
diff_x0 = xp0-x0

%---------------------------------
% simulate the models
%---------------------------------
Tt = 95;
Ti = 100;
T = Ti+Tt;
replications =200;
%Vectors of param values to generate matrix with shocks
sP1 = [rho_eps_tb sigma_tb_bar rho_sigma_tb eta_tb];
sP2 = [rho_eps_r sigma_r_bar rho_sigma_r eta_r r_bar 0];
[innovations,~] = my_innovations(T+1,replications,sP1,sP2);
% preallocate
R_simul_full = zeros(replications,Tt+1);

x_simul=zeros(model.n_x,T);
y_simul=zeros(model.n_y,T);
R_simul=zeros(model.n_y,T);
simul_data = zeros(3,T);
simul_moments = zeros(14,replications);

x_simul(:,1)=x0;

%option=1; % compute only simulated variables
option=2; % compute model residuals

% Simulation of exogenous states from perturbation (to make results
% comparable)
v_eps_r_tf = eps_r_t_1(:,2:end);
v_eps_tb_tf = eps_tb_t_1(:,2:end);

nexports = zeros(T,replications);

for n=1:replications
    
    x_simul(:,1)=x0;
    shocks = innovations(:,:,n);
    
    for t=1:T
        xt=x_simul(:,t);
        up0=shocks(:,t);
                
        % Option 2 - compute also model residuals
        [Rt,yt,~,auxt]=residual(coeffs,xt,params,c0,nodes,weights);
        
        y_simul(:,t)=yt;
        %keyboard;
        if t<T
        x_simul(:,t+1)=evalPhi(xt,yt,up0,params);
        x_simul(6,t+1)=v_eps_r_tf(n,t+1);
        x_simul(7,t+1)=v_eps_tb_tf(n,t+1);
                
        end
        R_simul(:,t)=Rt;
        lCt = auxt(1,1);lYt = auxt(7,1);lIft = auxt(13,1);
        simul_data(:,t) = [lCt;lYt;lIft];
               
    end
    
    consumption = exp(simul_data(1,:))';
    output = exp(simul_data(2,:))';
    investment = exp(simul_data(3,:))';
    %keyboard;
    nexports(:,n) = 100*(output - consumption - investment)./output;
    
    %keyboard;
    simul_data = simul_data(:,Ti:end);
    R_simul_full(n,:) = R_simul(1,Ti:end);
    simul_moments(:,n) = quarterly_moments_tp_full(simul_data);
        
end
%Exclude simulations with too large net exports (in abs. value)
cutoff = 95;
nexports = nexports(Ti:end,:);
max_nexports = max(abs(nexports),[],1);
v_netexports = abs(reshape(nexports,[1,replications*(Tt+1)]));
simul_moments = simul_moments(:,max_nexports<=prctile(v_netexports,cutoff));

std_y = simul_moments(1,:);
std_c_y = simul_moments(2,:);
std_i_y = simul_moments(3,:);
std_nx_y = simul_moments(4,:);
corr_nxy_y = simul_moments(7,:);

%average moments
m_std_y = mean(std_y);
m_std_c_y = mean(std_c_y);
m_std_i_y = mean(std_i_y);
m_std_nx_y = mean(std_nx_y);
m_corr_nxy_y = mean(corr_nxy_y);

R_simul_full1 = R_simul_full';
R_simul_full1 = R_simul_full1(:,max_nexports<=prctile(v_netexports,cutoff));

errors = log10(abs(R_simul_full1));

simulated_moments = zeros(8,1);
simulated_moments(1)= m_std_y*100;
simulated_moments(2)= m_std_c_y;
simulated_moments(3)= m_std_i_y;
simulated_moments(4)= nx_y_sss;
simulated_moments(5)= m_std_nx_y*100;
simulated_moments(6)= m_corr_nxy_y;
simulated_moments(7)= mean(mean(errors));
simulated_moments(8)= max(max(errors));

simulated_moments





