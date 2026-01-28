%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, solution around the DSS
% Computes deterministic debt transition (Figure 1, panel (a) in the
% paper), Euler equation errors from simulations, and ergodic distribution
% of external debt (Figure 1, panel (c) in the paper)
%---------------------------------------------------------------------------

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
xsss = mxt(:,Tsss+1:end);
ysss = myt(:,Tsss+1:end);

y_sss = exp(ysss(5));
c_sss = exp(ysss(1));
i_sss = exp(ysss(7));
D_sss = ysss(8);
r_sss = ysss(3);

%SSS net exports ratio
nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss
%SSS net exports ratio (alternative way to compute it)
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) )/y_sss
%Annualized external debt at DSS
nfa_dss = 100*nyss(8)/(12*exp(nyss(5)))
%Annualized external debt at SSS
nfa_sss = 100*D_sss/(12*y_sss)

nfa = 100*myt(8,:)./(12*exp(myt(5,:)));

%--------------------------------------------------------------------------
% Simulations to compute Euler errors and ergodic distribution of external
% debt
%--------------------------------------------------------------------------

%Prepare objects to compute Euler errors
P2 = [betta nu phi_d D_bar];
%Sigma = eye(n_e);
Sigma = M2f;
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_2(n_e,Sigma);

%Simulate and compute Euler errors
T0 = 1;
Tt = 95;
Ti = 100;
T = Ti+Tt;
replications =200;
errors =zeros(replications,T+1);
nfa_t = zeros(replications,T+1);


for t=1:replications
    shocks = innovations(:,:,t);
    %[yt,xt]=simul(xsss,shocks,nyss,nxss,eta,derivs,3,0,model);
    [yt,xt]=simul_mod_pruning3_dss(xsss,shocks,nyss,nxss,eta,derivs); %my pruning
    nfa_t(t,:) = 100*yt(8,:)./(12*exp(yt(5,:)));
    for s=1:T+1
        
        errors(t,s) =  log10(euler_errors_fgru_pac(P2,nxss,nyss,derivs,xt(:,s),eta,epsi_nodes,weight_nodes,approx1));
        
    end
    
end

nfa_t= nfa_t(:,Ti+1:end);
errors = errors(:,Ti+1:end);
my_mean_errors = (mean(mean(errors)))
my_max_errors = max(max(errors))

%Compute kernel density functions for NFA
%keyboard;
size_nfa = size(nfa_t);
nfa_t_res = reshape(nfa_t',size_nfa(1)*size_nfa(2),1);
pd_obj_b1 = fitdist(nfa_t_res,'Kernel','Kernel','epanechnikov');
XXb1=min(nfa_t_res):1:max(nfa_t_res);
pd_b_points1 = pdf(pd_obj_b1,XXb1);
figure;
subplot(1,2,1);
plot(nfa(1:600));xlabel('Time $t$');
subplot(1,2,2);
plot(XXb1,pd_b_points1,'b');xlabel('NFA');

