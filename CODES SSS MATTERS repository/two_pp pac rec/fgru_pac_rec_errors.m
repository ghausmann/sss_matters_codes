%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC Auxiliary model (SSS recalibrated "rec" solution)
% Computes Euler equation errors from simulations and the ergodic
% distribution of NFA
%---------------------------------------------------------------------------

%Parameters for aux. DSS (they're always fixed)
my_dss_params = [alppha r_bar delta etapar nu sigma_r_bar sigma_tb_bar];
% Make a vector of parameter values:
params=[r_bar, rho_eps_r, sigma_r_bar, rho_sigma_r, eta_r, rho_eps_tb, ...
    sigma_tb_bar, rho_sigma_tb, eta_tb, delta, alppha, nu, rho_x, betta, ...
    phipar, sigma_x, D_bar, rho_eps, etapar, phi_d, C_bar, psibetta, psi_k, psi_i];
% Make eta
% Define eta as in Schmitt-Grohe and Uribe (2004):
eta=[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;eta_r 0 0 0 0;0 eta_tb 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 exp(sigma_x);0 0 0 0 0];

%--------------------------------------------------------------------------
% Solving the model around the SSS 
%--------------------------------------------------------------------------

%Calibration of SSS (implied betta in model of interest, log(Ksss), excess
%depreciation in auxiliary model)
sss_sol = compute_betta_fgru_3p(D1,model,params,my_dss_params,M,eta,eps_ind,approx0,[-0.01 0.95*log(Kstar0) 0])
params(symparams=='psibetta') = sss_sol(1);
lK1 = sss_sol(2);
%Check that the SSS calibration delivers the targeted SSS
check_sss = compute_sss_fgru_3p(model,params,my_dss_params,M,eta,eps_ind,approx0,[D1 lK1 0])

%Compute aux. DSS and change parameter values ias needed
[nxss,nyss,psi_k,C_bar] = my_dss_aux([D1 sss_sol(2:3)],my_dss_params);
params(symparams=='D_bar') = D1;
params(symparams=='psi_k') = psi_k;
params(symparams=='C_bar') = C_bar;
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
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) + 0.5*phi_d*(D_sss-D1)^2   )/y_sss
%Annualized external debt at SSS
nfa_sss = 100*D_sss/(12*y_sss)

%--------------------------------------------------------------------------
% Simulations to compute Euler errors
%--------------------------------------------------------------------------

%Prepare objects to compute Euler errors
betta_moi = (1/(1+r_bar))*(((1+sss_sol(1))));
P2 = [betta_moi nu phi_d D_bar];
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

    [yt,xt]=simul_mod_pruning3(x0,shocks,nyss,nxss,eta,derivs);

    for s=1:T+1

        errors(t,s) =  log10(euler_errors_fgru_pac(P2,nxss,nyss,derivs,xt(:,s),eta,epsi_nodes,weight_nodes,approx1));
        
    end

    nfa_t(t,:) = 100*yt(8,:)./(12*exp(yt(5,:)));
end

errors = errors(:,Ti+1:end);
nfa_t= nfa_t(:,Ti+1:end);
my_mean_errors = (mean(mean(errors)))
my_max_errors = max(max(errors))

%Compute kernel density functions for NFA
size_nfa = size(nfa_t);
nfa_t_res = reshape(nfa_t',size_nfa(1)*size_nfa(2),1);
pd_obj_b1 = fitdist(nfa_t_res,'Kernel','Kernel','epanechnikov');
XXb1=min(nfa_t_res):1:max(nfa_t_res);
pd_b_points1 = pdf(pd_obj_b1,XXb1);
figure;plot(XXb1,pd_b_points1,'b');