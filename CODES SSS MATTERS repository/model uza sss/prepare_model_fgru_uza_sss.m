%
% FGRU Small open-economy model
% UZAWA version of the model
% It works both for standard perturbation around DSS, and two-parameter
% perturbation around SSS.
% It uses Levintal's perturbation toolbox
% INSTRUCTIONS: Run this code just once
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

%Note: ug_r and ug_tb are auxiliary states so that we can introduce
%uncertainty shocks with the canonical form (see right below):
%ug_r_{t+1}=u_r_{t+1}, E(u_r_{t+1})=0.
%Note: eps is the second perturbation object with law of motion
% eps_{t+1}=rho_eps*eps_{t}

clear;

%Start fresh with the search path, and activate the required subfolders
my_root = locate_root(mfilename('fullpath'));
remove_subfolders(my_root)
methods_path = genpath(fullfile(my_root,'sss_matters_methods'));
addpath(methods_path);

% Step 1: Define the model by symbolic variables.

% Declare parameters
syms r_bar rho_eps_r sigma_r_bar rho_sigma_r eta_r rho_eps_tb ...
    sigma_tb_bar rho_sigma_tb eta_tb delta alppha nu rho_x betta ...
    phipar sigma_x D_bar etapar rho_eps C_bar phi_d ...
    psibetta psi_i psi_k real;
% Declare state variables:
% today
syms sigma_r sigma_tb eps_r eps_tb X D K eps I ug_tb ug_r real;
% tomorrow
syms sigma_rp sigma_tbp eps_rp eps_tbp Xp Dp Kp epsp Ip ug_tbp ug_rp real;
% Declare control variables: 
% today
syms C H r lambda Y phi If Df real;
% tomorrow
syms Cp Hp rp lambdap Yp phip Ifp Dfp real;

%Collect state variables: 
x=[D, K, I, sigma_r, sigma_tb, eps_r, eps_tb, ug_r, ug_tb, X, eps]; % today 
xp=[Dp, Kp, Ip, sigma_rp, sigma_tbp, eps_rp, eps_tbp, ug_rp, ug_tbp, Xp, epsp]; % tomorrow

%Control variable: log consumption
y=[C, H, r, lambda, Y, phi, If, Df]; % today
yp=[Cp, Hp, rp, lambdap, Yp, phip, Ifp, Dfp]; % tomorrow


%Collect parameters
symparams=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, rho_eps, etapar, phi_d, C_bar, psibetta, psi_k, psi_i];
%Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.

%Note: here logz is a pre-determined state, just like bonds
Phi = [ (1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r; (1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb; ...
    rho_eps_r*eps_r+exp(sigma_r)*ug_r; rho_eps_tb*eps_tb+exp(sigma_tb)*ug_tb;  ...
    0; 0; rho_x*X;rho_eps*eps];

% Define eta as in Schmitt-Grohe and Uribe (2004):
eta=[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;eta_r 0 0 0 0;0 eta_tb 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 exp(sigma_x);0 0 0 0 0];

%budget constraint
f1 = exp(Y) - exp(K)^(alppha)*(exp(X)*exp(H))^(1-(alppha)); %Eq 2.1 in paper

f2 = exp(Kp) - (1-(delta+psi_i*(1-eps^2)))*exp(K) - (1-phipar/2*(exp(If)/exp(I)-1)^2)*exp(Ip); %Eq 2.25 in paper

f3 = exp(Y)-exp(C)-exp(If) - (D) + (Dp)/(1+r); %Eq 2.18 in paper

f4 = r - ( r_bar  ) - (rho_eps_tb*eps_tb+exp(sigma_tb)*ug_tb) - (rho_eps_r*eps_r+exp(sigma_r)*ug_r); %Eq 2.4 in paper (combined with 2.10 and 2.11)

f5 = (exp(C))^(-nu) - exp(lambda); %Eq 2.5 in paper

f6 = exp(lambda)*( 1/(1+ r) ) -  ...
  betta*(1+psibetta*eps^2)*((exp(C)/C_bar)^-phi_d)*exp(lambdap); %Eq 2.15 in paper

f7 = -exp(phi)+betta*(1+psibetta*eps^2)*((exp(C)/C_bar)^-phi_d)*...
    ((1-(delta+psi_k*(1-eps^2)))*exp(phip)+(alppha)*exp(Yp)/exp(Kp)*exp(lambdap)); %Capital Euler eq. in Appendix A.

f8 = (C_bar^-nu)*(exp(C)^nu)*(exp(H)^etapar) - (1-(alppha))*exp(Y); %Eq 2.8 in paper

f9 = exp(phi)*(1-phipar/2*((exp(If)-exp(I))/exp(I))^2-phipar*exp(If)/exp(I)*((exp(If)-exp(I))/exp(I))) ...
    + betta*(1+psibetta*eps^2)*((exp(C)/C_bar)^-phi_d)*exp(phip)*phipar*(exp(Ifp)/exp(If))^2*((exp(Ifp)-exp(If))/exp(If)) -exp(lambda); %Eq 2.17 in paper

f10 = exp(Ip) - exp(If); %Definition of If

f11 = Df - Dp; %Definition of Df

f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11]; % all conditions

% Choose approximation order
approx=3; % fifth order solution.

% call differentiate_dsge

model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
%save('model')
save('files\model.mat')