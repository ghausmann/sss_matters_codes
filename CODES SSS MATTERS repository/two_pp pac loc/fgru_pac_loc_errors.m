%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Auxiliary model (SSS fixed), approximation around the SSS
% Computes Euler equation errors from simulations
%---------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Solving the model around the SSS 
%--------------------------------------------------------------------------

%Solve for the SSS
sss_sol = compute_sss_fgru_3_locp_num(model,params,my_dss_params,M,eta,eps_ind,approx0,[1 log(Kstar0) 0])
%Recompute DSS of auxiliary model
[nxss,nyss,psi_k] = my_num_dss(sss_sol,my_dss_params);

params(symparams=='Da') = sss_sol(1);
params(symparams=='psi_k') = psi_k;
params(symparams=='psi_i') = sss_sol(3);

% Compute the perturbation solution:
%DSS check
mycheck0 = double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
mycheck=sum(abs(mycheck0))
%Calculate derivatives of the policy functions abouts the DSS
derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx1,algo);

%------------------------------------------------------------------------------
% SSS results
%------------------------------------------------------------------------------

x0=nxss; % start at the steady state
derivs.hx(eps_ind,eps_ind) = 1;
x0(eps_ind) = 1; %impose the model of interest (epsilon=1)

%Evaluate decision rules for states and controls
mh_sss = dr_ht(derivs,nxss,approx0,(x0-nxss));
mg_sss = dr_gt(derivs,nyss,approx0,(x0-nxss));
y_sss = exp(mg_sss(y=='Y'));
c_sss = exp(mg_sss(y=='C'));
i_sss = exp(mg_sss(y=='If'));
D_sss = mg_sss(y=='Df');
r_sss = mg_sss(y=='r');

%SSS net exports ratio
nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss
%SSS net exports ratio (alternative way to compute it)
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) + 0.5*phi_d*(D_sss-D_bar)^2   )/y_sss
%Punishment for not holding a D1 exogenous debt level
crazy_punish = 100*(  0.5*phi_d*(D_sss-D_bar)^2   )/y_sss
%Annualized external debt at DSS
nfa_dss = 100*nyss(8)/(12*exp(nyss(5)))
%Annualized external debt at SSS
nfa_sss = 100*D_sss/(12*y_sss)

%--------------------------------------------------------------------------
% Simulations to compute Euler errors
%--------------------------------------------------------------------------

%Prepare objects to compute Euler errors
betta_moi = (1/(1+r_bar));
P2 = [betta_moi nu phi_d D_bar];
Sigma = M2f;
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_2(5,Sigma);

%Simulate and compute Euler errors
T0 = 1;
Tt = 95;
Ti = 100;
T = Ti+Tt;
replications =200;
errors =zeros(replications,T+1);
output = zeros(replications,T+1);
consumption = zeros(replications,T+1);
investment = zeros(replications,T+1);
nexports = zeros(replications,T+1);
nfa_t = zeros(replications,T+1);

for t=1:replications
    
    shocks = innovations(:,:,t);
    if act_pruning==1
        [yt,xt]=simul_mod_pruning3(x0,shocks,nyss,nxss,eta,derivs);
    else
        [yt,xt]=simul(x0,shocks,nyss,nxss,eta,derivs,approx1,0,model);
    end
    
    for s=1:T+1
        
        errors(t,s) =  log10(euler_errors_fgru_pac(P2,nxss,nyss,derivs,xt(:,s),eta,epsi_nodes,weight_nodes,approx1));
        
    end
    
    output(t,:) = exp(yt(5,:));
    consumption(t,:) = exp(yt(1,:));
    investment(t,:) = exp(yt(7,:));
    nexports(t,:) = 100*(output(t,:) - consumption(t,:) - investment(t,:))./output(t,:);
    nfa_t(t,:) = 100*yt(8,:)./(12*exp(yt(5,:)));
end

errors = errors(:,Ti+1:end);
nfa_t= nfa_t(:,Ti+1:end);
nexports= nexports(:,Ti+1:end);

%If the percentile cutoff is below 100, this will exclude simulations with
%large absolute export values
max_nexports = max(abs(nexports),[],2);
v_netexports = abs(reshape(nexports,[1,replications*(Tt+1)]));
errors_win = errors(max_nexports<=prctile(v_netexports,pctile),:);
my_mean_errors = (mean(mean(errors_win)))
my_max_errors = max(max(errors_win))
 