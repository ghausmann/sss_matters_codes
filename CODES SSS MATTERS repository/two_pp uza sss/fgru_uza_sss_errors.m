%---------------------------------------------------------------------------
% FGRU Small open-economy model
% Uzawa model, approximation around the SSS
% Computes Euler equation errors from simulations
%---------------------------------------------------------------------------

%Make vector of parameter values and matrices of innovations
routine_sss_tasks;

%--------------------------------------------------------------------------
% Solving the model around the SSS 
%--------------------------------------------------------------------------

%Calibration of SSS (implied betta in model of interest, log(Ksss), excess
%depreciation in auxiliary model)
sss_sol = compute_betta_fgru_3p(D1,model,params,my_dss_params,M,eta,eps_ind,approx0,[-0.01 0.9*log(Kstar0) 0])
params(symparams=='psibetta') = sss_sol(1);
lK1 = sss_sol(2);
%Check that the SSS calibration delivers the targeted SSS
check_sss = compute_sss_fgru_3p(model,params,my_dss_params,M,eta,eps_ind,approx0,[D1 lK1 0])

%Compute aux. DSS and change parameter values as needed
[nxss,nyss,psi_k,C_bar] = my_dss_aux([D1 sss_sol(2:3)],my_dss_params)
params(symparams=='D_bar') = D1;
params(symparams=='psi_k') = psi_k;
params(symparams=='C_bar') = C_bar;
params(symparams=='psi_i') = psi_i;
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

nx_y_sss = 100*(y_sss - c_sss - i_sss)/y_sss
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss))  )/y_sss
nfa_sss = 100*D_sss/(12*y_sss)

%--------------------------------------------------------------------------
% Simulations to compute Euler errors
%--------------------------------------------------------------------------

%Prepare objects to compute Euler errors
betta_moi = (1/(1+r_bar))*(((1+sss_sol(1))));
P2 = [betta_moi nu phi_d C_bar];
Sigma = M2f;
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_2(n_e,Sigma);

%Simulate and compute Euler errors
T0 = 1;
Tt = 95;
Ti = 100;
T = Ti+Tt;
replications =200;
errors =zeros(replications,T+1);

for t=1:replications
    
    shocks = innovations(:,:,t);
    if act_pruning==1
        [yt,xt]=simul_mod_pruning3(x0,shocks,nyss,nxss,eta,derivs);
    else
        [yt,xt]=simul(x0,shocks,nyss,nxss,eta,derivs,approx1,0,model);
    end
    
    for s=1:T+1
        
        errors(t,s) =  log10(euler_errors_fgru_uzawa(P2,nxss,nyss,derivs,xt(:,s),eta,epsi_nodes,weight_nodes,approx1));
        
    end
        
end

errors = errors(:,Ti+1:end);
my_mean_errors = (mean(mean(errors)))
my_max_errors = max(max(errors))
