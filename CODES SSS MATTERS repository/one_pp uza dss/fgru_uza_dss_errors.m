%---------------------------------------------------------------------------
% FGRU Small open-economy model
% UZAWA model, solution around DSS
% Computes deterministic debt transition and Euler equation errors 
% from simulations(Table 5 in paper)
%---------------------------------------------------------------------------

D_bar = D1;
betta = 1/(1+(r_bar));   %discount factor

%Compute the perturbation solution
routine_sol_dss_tasks_uza;

%------------------------------------------------------------------------------
% simulate the model for Tsss periods (FIRST RUN to converge to an implied SSS)
%------------------------------------------------------------------------------
x0=nxss; % start at the steady state
Tsss = 2000;
shocks0 = zeros(n_e,Tsss); % a bunch of zeros
[myt,mxt]=simul(x0,shocks0,nyss,nxss,eta,derivs,3,0,model);
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
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss))  )/y_sss
%Annualized external debt at DSS
nfa_dss = 100*nyss(8)/(12*exp(nyss(5)))
%Annualized external debt at SSS
nfa_sss = 100*D_sss/(12*y_sss)

nfa = 100*myt(8,:)./(12*exp(myt(5,:)));
figure;plot(nfa);xlabel('time');title('convergence NFA over time');

%--------------------------------------------------------------------------
% Simulations to compute Euler errors
%--------------------------------------------------------------------------

%Prepare objects to compute Euler errors
P2 = [betta nu phi_d C_bar];
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
[innovations,rates] = my_innovations(T,replications,sP1,sP2);
nfa_t = zeros(replications,T+1);

for t=1:replications
    shocks = innovations(:,:,t);
    %[yt,xt]=simul(xsss,shocks,nyss,nxss,eta,derivs,3,0,model);
    [yt,xt]=simul_mod_pruning3_dss(xsss,shocks,nyss,nxss,eta,derivs); %my pruning
    nfa_t(t,:) = 100*yt(8,:)./(12*exp(yt(5,:)));
    for s=1:T+1
        
        errors(t,s) =  log10(euler_errors_fgru_uzawa(P2,nxss,nyss,derivs,xt(:,s),eta,epsi_nodes,weight_nodes,approx1));
        
    end
    
end
    
errors = errors(:,Ti+1:end);
nfa_t= nfa_t(:,Ti+1:end);
my_mean_errors = (mean(mean(errors)))
my_max_errors = max(max(errors))
