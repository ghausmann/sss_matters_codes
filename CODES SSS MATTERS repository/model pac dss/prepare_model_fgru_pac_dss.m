%
% FGRU Small open-economy model
% Basic version as in AER 2011 paper
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
% eps_{t+1}=rho_eps*eps_{t} (but here it does not play any role).

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
    phipar sigma_x D_bar etapar rho_eps C_bar phi_d real;
% Declare state variables:
% today (note: D and I are debt and investment at time t, treated as states)

syms sigma_r sigma_tb eps_r eps_tb X D K eps I ug_tb ug_r real;
% tomorrow
syms sigma_rp sigma_tbp eps_rp eps_tbp Xp Dp Kp epsp Ip ug_tbp ug_rp real;
% Declare control variables: 
% today (note: Df and If are debt and investment at time t, treated as controls)
syms C H r lambda Y varphi If Df real;
% tomorrow
syms Cp Hp rp lambdap Yp varphip Ifp Dfp real;

%Collect state variables: 
x=[D, K, I, sigma_r, sigma_tb, eps_r, eps_tb, ug_r, ug_tb, X, eps]; % today 
xp=[Dp, Kp, Ip, sigma_rp, sigma_tbp, eps_rp, eps_tbp, ug_rp, ug_tbp, Xp, epsp]; % tomorrow

%Collect control variables:
y=[C, H, r, lambda, Y, varphi, If, Df]; % today
yp=[Cp, Hp, rp, lambdap, Yp, varphip, Ifp, Dfp]; % tomorrow

%Collect parameters
symparams=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, rho_eps, etapar, phi_d, C_bar];

%Define the function Phi, which is the lower block of h(x).
%This block is the expected value of the exogenous state variables.
%Includes equations (2.13), (2.12), (2.11), (2.10), and (2.14) in paper 
Phi = [ (1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r; (1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb; ...
    rho_eps_r*eps_r+exp(sigma_r)*ug_r; rho_eps_tb*eps_tb+exp(sigma_tb)*ug_tb;  ...
    0; 0; rho_x*X;rho_eps*eps];
%Define eta as in Schmitt-Grohe and Uribe (2004):
eta=[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;eta_r 0 0 0 0;0 eta_tb 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 exp(sigma_x);0 0 0 0 0];

%Equilibrium conditions
f1 = exp(Y) - exp(K)^(alppha)*(exp(X)*exp(H))^(1-(alppha)); %Eq 2.1 in paper

f2 = exp(Kp) - (1-delta)*exp(K) - (1-phipar/2*(exp(If)/exp(I)-1)^2)*exp(Ip); %Eq 2.2 in paper

f3 = exp(Y)-exp(C)-exp(If) - (D) + (Dp)/(1+r) - 0.5*phi_d*(Dp-D_bar)^2; %Eq 2.3 in paper

f4 = r - ( r_bar ) - (rho_eps_tb*eps_tb+exp(sigma_tb)*ug_tb) - (rho_eps_r*eps_r+exp(sigma_r)*ug_r); %Eq 2.4 in paper (combined with 2.10 and 2.11)

f5 = (exp(C))^(-nu) - exp(lambda); 

f6 = exp(lambda)*( 1/(1+ r) - phi_d*(Dp-D_bar)    ) -  ...
  betta*exp(lambdap); %Eq 2.6 in paper

f7 = -exp(varphi)+betta*((1-delta)*exp(varphip)+(alppha)*exp(Yp)/exp(Kp)*exp(lambdap)); %Eq 2.7 in paper

f8 = (C_bar^-nu)*(exp(C)^nu)*(exp(H)^etapar) - (1-(alppha))*exp(Y); %Eq 2.8 in paper

f9 = exp(varphi)*(1-phipar/2*((exp(If)-exp(I))/exp(I))^2-phipar*exp(If)/exp(I)*((exp(If)-exp(I))/exp(I))) ...
    + betta*exp(varphip)*phipar*(exp(Ifp)/exp(If))^2*((exp(Ifp)-exp(If))/exp(If)) -exp(lambda); %Eq 2.9 in paper

f10 = exp(Ip) - exp(If); %Definition of If

f11 = Df - Dp; %Definition of Df 

f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11]; % all conditions

% Choose approximation order

approx=5; % fifth order solution.

% call differentiate_dsge

model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
save('files\model.mat')