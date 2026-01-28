/*
 *  FGRU Small open-economy model
 *  
 */

//******************************Dynare definitions************************************************

%Make sure declaration order matches DR order
%Variables Dp and Ip are just the leads of debt and investment (at time t,Dp=D and Ip=I)
%Note: eps is the second perturbation object with law of motion
%eps_{t+1}=rho_eps*eps_{t}, but here it does not play any role (does not enter in any equilibrium condition).
var C H r Dp sigma_r sigma_tb eps_r eps_tb X D K I eps lambda Y phi Ip;
%Order of innovations same as in codes using Levintal's toolbox
varexo u_sigma_r u_sigma_tb u_r u_tb u_x;

parameters r_bar rho_eps_r sigma_r_bar rho_sigma_r eta_r rho_eps_tb 
           sigma_tb_bar rho_sigma_tb eta_tb delta alppha nu rho_x betta
           phipar sigma_x D_bar rho_eps etapar phi_d C_bar;
            
//******************************Dynare parameter initialization************************************************

// Common parameters:
rho_eps_tb=0.95;
sigma_tb_bar=-8.06;     
rho_sigma_tb=0.94;
eta_tb=0.13;
nu=5;                   
etapar = 1001;
delta=0.014;            
alppha = 0.32;          
rho_x=0.95;
rho_eps = 1; %persistence of new perturbation state             
            
// Country-specific parameters (Argentina):
rho_eps_r=0.97;
sigma_r_bar=-5.71;
rho_sigma_r=0.94;
eta_r=0.46;
r_bar=(0.02);
betta = 1/(1+(r_bar));  %discount factor.

// Estimated parameters
phi_d = 7.4686e-05; %PAC parameter 
phipar=64.8784; %capital adjustment costs
sigma_x=-3.2567; %std. of innovations to productivity
D_bar = 25.7968; %exogenous debt level

// DSS in the auxiliary model
Hstar = ( (1-alppha)*(alppha/((r_bar)+delta))^(alppha/(1-alppha)))^(1/(etapar-1));
Kstar = Hstar*(alppha/((r_bar)+delta))^(1/(1-alppha));
Istar = delta*Kstar;
Ystar = Kstar^alppha*Hstar^(1-alppha);
Cstar = D_bar/(1+(r_bar))-D_bar+Ystar-Istar;
C_bar = Cstar;
lambdastar = (Cstar)^(-nu);
phistar = lambdastar;

//******************************Dynare model equations************************************************

model;

%Eq. 2.1 in paper:
exp(Y)=exp(K(-1))^(alppha)*(exp(X)*exp(H))^(1-(alppha)); 
%Eq. 2.2 in paper:
exp(K) = (1-(delta))*exp(K(-1)) + (1-phipar/2*(exp(I)/exp(I(-1))-1)^2)*exp(I);
%Eq. 2.3 in paper:
exp(Y)-exp(C)-exp(I) = D(-1) - D/(1+r) + 0.5*phi_d*(D-D_bar)^2;
%Eq. 2.4 in paper:
r=(r_bar) + eps_tb + eps_r;

%Eq. 2.5 in paper:
(exp(C))^(-nu)=exp(lambda);
%Eq. 2.6 in paper:
exp(lambda)*(1/(1+ r) - phi_d*(D-D_bar) ) = betta*exp(lambda(+1));
%Eq. 2.7 in paper:
exp(phi) = betta*( ( 1-( delta ) )*exp(phi(+1))+(alppha)*(exp(Y(+1))/exp(K))*exp(lambda(+1)));
%Eq. 2.8 in paper:
(C_bar^-nu)*(exp(C)^nu)*(exp(H)^etapar)=(1-(alppha))*exp(Y);
%Eq. 2.9 in paper:
exp(phi)*(1-phipar/2*((exp(I)-exp(I(-1)))/exp(I(-1)))^2-phipar*exp(I)/exp(I(-1))*((exp(I)-exp(I(-1)))/exp(I(-1))))
+ betta*exp(phi(+1))*phipar*(exp(Ip(+1))/exp(Ip))^2*((exp(Ip(+1))-exp(Ip))/exp(Ip))=exp(lambda);

%Definitions of leads Dp and Ip:
Dp = D;
Ip = I; 

//Exogenous laws of motion
%Eqs. 2.10 and 2.11 in paper
eps_tb=rho_eps_tb*eps_tb(-1)+exp(sigma_tb)*u_tb;   
eps_r=rho_eps_r*eps_r(-1)+exp(sigma_r)*u_r;   
%Eqs. 2.12 and 2.13 in paper
sigma_tb=(1-rho_sigma_tb)*sigma_tb_bar+rho_sigma_tb*sigma_tb(-1)+eta_tb*u_sigma_tb;
sigma_r=(1-rho_sigma_r)*sigma_r_bar+rho_sigma_r*sigma_r(-1)+eta_r*u_sigma_r;
%Eq. 2.14 in paper
X=rho_x*X(-1)+exp(sigma_x)*u_x;
%Law of motion for new perturbation state:
eps = rho_eps*eps(-1);

end;


//******************************Initial values for steady state computation******************************************
initval;
%Exogenous states
sigma_tb=sigma_tb_bar;
sigma_r=sigma_r_bar;
eps_r=0;
eps_tb=0;
X = 0;
eps = 0;

%innovations
u_x = 0;
u_r = 0;
u_tb = 0;
u_sigma_tb = 0;
u_sigma_r = 0;

% steady states
C = log(Cstar);
K = log(Kstar);
lambda = log(lambdastar);
H = log(Hstar);
phi = log(phistar);
D=D_bar;
I = log(Istar);
Y = log(Ystar);
r=r_bar;
Dp = D_bar;
Ip = log(Istar);
end;

//******************************Set covariance matrix of structural shocks to identity**************************************
shocks;
var u_x; stderr 1;
var u_r; stderr 1;
var u_tb; stderr 1;
var u_sigma_tb; stderr 1;
var u_sigma_r; stderr 1;
end;

//******************************Compute dynamic solution**************************************

steady;
check;
//get decision rules at order 3;
stoch_simul(order=3,periods=0,irf=0); 

