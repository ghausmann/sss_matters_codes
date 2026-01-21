% This function implements the SSS algorithm as described in Appendix C of
% the paper, using the auxiliary model of Section 2.2 in the paper, and set
% in "calibration" mode (by keeping the exogenous debt level fixed, and
% finding the discount factor consistent with that target).

function y = calib_sss_fgru_3p(x,D1,model,params,my_dss_params,M,eta,eps_ind,approx)

load('model','symparams');

algo='vectorize'; % Simple vectorization.
%algo='dlyap'; % Hessenberg-Schur algorithm.
%algo='gensylv'; % Kamenik algorithm

lK1 = x(2); % approximation point for capital
psi_i = x(3);
delta = params(symparams=='delta');

%Recompute DSS of auxiliary model
[nxss,nyss,psi_k,C_bar] = my_dss_aux([D1 x(2:3)],my_dss_params);

params(symparams=='psibetta') = x(1);
params(symparams=='D_bar') = D1;
params(symparams=='psi_k') = psi_k;
params(symparams=='C_bar') = C_bar;
params(symparams=='psi_i') = psi_i;

%Compute up to k-order derivatives
%use Levintal's function solve_dsge.m
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

%Evaluate at DSS and sigma=epsilon=1
x0=nxss; % start at the steady state
if eps_ind>0
    x0(eps_ind) = 1; % impose the model of interest
    derivs1.hx(eps_ind,eps_ind) = 1;
end

%Compute the Taylor-series policy rule for bonds
%uses my own function dr_ht
x1=dr_ht(derivs1,nxss,approx,(x0-nxss));

%Return h1(DSS,1,1)-b_bar 
eq1 = x1(1) - D1;
eq2 = x1(2) - (lK1);
eq3 = x1(3) - log((delta+psi_i)*exp(lK1));
%keyboard;
%eq3 = x1(3) - log(I1);

y = [eq1;eq2;eq3];