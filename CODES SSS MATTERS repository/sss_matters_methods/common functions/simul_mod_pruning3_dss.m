function [yt,xt]=simul_mod_pruning3_dss(x0,shocks,nyss,nxss,eta,derivs)
% [yt,xt]=simul_mod_pruning3_dss(x0,shocks,nyss,nxss,eta,derivs) is a
% heavily modified version of the original function simul.m by Oren
% Levintal (2016). You can use it to simulate a two-parameter perturbation
% model where third-order policy functions are approximated around the
% deterministic steady-state of the state variables.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

lx = numel(nxss);
x0(lx) = 0;
x_0 = nxss;

% Compute jacobians accounting for third-order perturbation effects
derivs_hx_c = my_jac_hx(derivs,nxss,3,0.0001);
derivs_gx_c = my_jac_gx(derivs,nxss,nyss,3,0.0001);

% Evaluate policy rules at (sigma,epsilon)=(1,1) to find SSS
nnyss1 = dr_gt(derivs,nyss,2,(x_0-nxss));
nnxss1 = dr_ht(derivs,nxss,2,(x_0-nxss));

T=size(shocks,2);
n_y=size(derivs.gx,1);
n_x=size(derivs.hx,1);
n_e=size(shocks,1);
shocks=[zeros(n_e,1),shocks,zeros(n_e,1)];

xt_f=zeros(n_x+1,T+2);
yt=zeros(n_y,T+2);
xt_s=zeros(n_x+1,T+2);
xt_rd=zeros(n_x+1,T+2);
xt_f(1:end-1,1)=x0-nxss;


for t=1:T+1
    x_f=xt_f(:,t);
    xt_f(1:end-1,t+1)= (nnxss1-nxss) + derivs_hx_c*x_f+eta*shocks(:,t+1);
    
    x_s=xt_s(:,t);
    x_f2=kron(x_f,x_f);
    xt_s(1:end-1,t+1)=  derivs_hx_c*x_s+derivs.hxx*x_f2/2;
    
    x_rd=xt_rd(:,t);
    x_f3=kron(x_f2,x_f);
    x_f_x_s=kron(x_f,x_s);
    xt_rd(1:end-1,t+1)= derivs_hx_c*x_rd+derivs.hxx*(2*x_f_x_s)/2+derivs.hxxx*x_f3/6;
    
    
    yt(:,t)=  (nnyss1-nyss) +  derivs_gx_c*(x_f+x_s+x_rd)+derivs.gxx*(x_f2+2*x_f_x_s)/2 ...
        +derivs.gxxx*(x_f3)/6;
    
end

yt=yt(:,1:T+1);
xt=xt_f(:,1:T+1) + xt_s(:,1:T+1) + xt_rd(:,1:T+1);


yt=yt+repmat(nyss,1,T+1);
xt=xt(1:end-1,:)+repmat(nxss,1,T+1);


