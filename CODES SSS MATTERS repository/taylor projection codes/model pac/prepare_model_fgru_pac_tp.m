%
% FGRU Small open-economy model
%

% The model is defined by:
%
% Ef(yp,y,xp,x)=0
% y=g(x);
% xp=h(x)+eta*ep
% h(x)=[tilh(x);Phi(x)];
%
% The functions f and Phi are known.
% The functions g and tilh are unknown.

% The equilibrium conditions include homotopy "tricks" to help the
% nonlinear solver to achieve convergence.

clear;

my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
methods_path = genpath(fullfile(my_root,'taylor projection codes\solution_methods'));
addpath(methods_path);

% Step 1: Define the model by symbolic variables.

% Declare parameters
syms r_bar rho_eps_r sigma_r_bar rho_sigma_r eta_r rho_eps_tb ...
    sigma_tb_bar rho_sigma_tb eta_tb delta alppha nu rho_x betta ...
    phipar sigma_x D_bar etapar C_bar phi_d risk_factor real;
% Declare state variables:
% today
syms sigma_r sigma_tb eps_r eps_tb X D K I u_tb u_r real;
% tomorrow
syms sigma_rp sigma_tbp eps_rp eps_tbp Xp Dp Kp Ip u_tbp u_rp real;
% Declare control variables: 
% today
syms C H r phi Df real;
% tomorrow
syms Cp Hp rp phip Dfp real;
%innovations;
syms up_sigma_r up_sigma_tb up_r up_tb up_x real;

%-------------------------------
% Define the auxiliary functions
%-------------------------------

syms If lC lH llambda lY lphi lI lIf lK lX IP real;
syms Ifp lCp lHp llambdap lYp lphip lIp lIfp lKp lXp IPp real;

Df_ = Dp;
lC_ = exp(C); lCp_ = exp(Cp);
lH_ = exp(H); lHp_ = exp(Hp);
llambda_ = lC^-nu; llambdap_ = lCp^-nu;
lphi_ = exp(phi); lphip_ = exp(phip);
lI_ = exp(I); lIp_ = exp(Ip);
lIf_ = exp(If); lIfp_ = exp(Ifp);
lK_ = exp(K); lKp_ = exp(Kp);
lX_ = exp(X); lXp_ = exp(Xp);
lY_ = (lK^(alppha)*(lX*lH)^(1-(alppha))); lYp_ = (lKp^(alppha)*(lXp*lHp)^(1-(alppha)));
IP_ = (lIf-lI)/lI; IPp_ = (lIfp-lIp)/lIp;

%-----------------------------------------------------------
% Vectors of auxiliary functions and corresponding variables
%-----------------------------------------------------------

auxfuns=[lC_;lCp_;lH_;lHp_;llambda_;llambdap_;lY_;lYp_;lphi_;lphip_;lI_;lIp_;lIf_;lIfp_;lK_;lKp_;lX_;lXp_;IP_;IPp_];
auxvars=[lC;lCp;lH;lHp;llambda;llambdap;lY;lYp;lphi;lphip;lI;lIp;lIf;lIfp;lK;lKp;lX;lXp;IP;IPp];

%Collect state variables: 
x=[D, K, I, sigma_r, sigma_tb, eps_r, eps_tb, X]; % today 
xp=[Dp, Kp, Ip, sigma_rp, sigma_tbp, eps_rp, eps_tbp, Xp]; % tomorrow

%Collect control variables:
y=[C, H, r, phi, If, Df]; % today
yp=[Cp, Hp, rp, phip, Ifp, Dfp]; % tomorrow

%-----------------
% Vector of shocks
%-----------------
shocks=[up_sigma_r,up_sigma_tb,up_r,up_tb,up_x];

%Collect parameters
symparams=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, etapar, phi_d, C_bar, risk_factor];

%Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.
Phi_fun = [Df;log((1-(delta))*(lK) + (1-phipar/2*(IP)^2)*lIf);If;
    (1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r + risk_factor*eta_r*up_sigma_r;
    (1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb + risk_factor*eta_tb*up_sigma_tb;
    rho_eps_r*eps_r+risk_factor*exp((1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r + risk_factor*eta_r*up_sigma_r)*up_r;
    rho_eps_tb*eps_tb+risk_factor*exp((1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb + risk_factor*eta_tb*up_sigma_tb)*up_tb;
    rho_x*X + risk_factor*exp(sigma_x)*up_x];
    

f1 = ( lC + lIf + (D) - (Dp)/(1+r) + 0.5*phi_d*(Dp-D_bar)^2 )/lY - 1; % Eq. 2.3 in paper
%
f2 = r - (r_bar) - eps_tb - eps_r; % Eq. 2.4 in paper
%
f3 = (( 1/(1+ r) - phi_d*(Dp-D_bar) )^-1)*betta*(llambdap/llambda) - 1; % Eq. 2.6 in paper
%
f4 = betta*((1-(delta))*lphip+(alppha)*lYp/lKp*llambdap)/lphi - 1; % Eq. 2.7 in paper
%
f5 = ((C_bar^-nu)*(lC^nu)*(lH^etapar))/((1-(alppha))*lY) - 1; % Eq. 2.8 in paper
%
f6 = lphi*(1-phipar/2*(IP)^2-phipar*(IP+1)*IP )/llambda ...
    + betta*lphip*phipar*(IPp+1)^2*( IPp )/llambda - 1; % Eq. 2.9 in paper

f=[f1;f2;f3;f4;f5;f6]; % all conditions

% Choose approximation order
order=3; % fifth order solution.

%----------------
% Call prepare_tp
%----------------

model=prepare_tp(f,Phi_fun,yp,y,xp,x,shocks,symparams,order,auxfuns,auxvars);

%-----------
% Save model
%-----------

save('model') % you will need this later
