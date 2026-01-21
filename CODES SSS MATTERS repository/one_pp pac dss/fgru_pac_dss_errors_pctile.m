%---------------------------------------------------------------------------
% FGRU Small open-economy model
% PAC model, solution around the DSS
% Computes Euler equation errors from simulations (to replicate numbers from Table 3 in the paper)
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
nx_y_sss_alt = 100*( D_sss*(1 - 1/(1+r_sss)) + 0.5*phi_d*(D_sss-D_bar)^2   )/y_sss
%Punishment for not holding a D1 exogenous debt level
crazy_punish = 100*(  0.5*phi_d*(D_sss-D_bar)^2   )/y_sss
%Annualized external debt at DSS
nfa_dss = 100*nyss(8)/(12*exp(nyss(5)))
%Annualized external debt at SSS
nfa_sss = 100*D_sss/(12*y_sss)

nfa = 100*myt(8,:)./(12*exp(myt(5,:)));

%--------------------------------------------------------------------------
% Simulations to compute Euler errors
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

%load('my_dss_error_5_data.mat');

errors =zeros(replications,T+1);
r_ap = zeros(replications,T+1);
nfa_t = zeros(replications,T+1);
nexports = zeros(replications,T+1);


for t=1:replications
    shocks = innovations(:,:,t);
    [yt,xt]=simul(xsss,shocks,nyss,nxss,eta,derivs,approx1,0,model); %no pruning
    nfa_t(t,:) = 100*yt(8,:)./(12*exp(yt(5,:)));

    output = exp(yt(5,:));
    consumption = exp(yt(1,:));
    investment = exp(yt(7,:));
    nexports(t,:) = 100*(output - consumption - investment)./output;


    for s=1:T+1
        
        errors(t,s) =  log10(euler_errors_fgru_pac(P2,nxss,nyss,derivs,xt(:,s),eta,epsi_nodes,weight_nodes,approx1));
        
    end
   % count = count+1
end

%Exclude simulations with large abs. values of net exports, beyond the
%pctile cutoff
cutoff = 95;
errors = errors(:,Ti+1:end);
nfa_t= nfa_t(:,Ti+1:end);
nexports= nexports(:,Ti+1:end);

max_nexports = max(abs(nexports),[],2);
v_netexports = abs(reshape(nexports,[1,replications*(Tt+1)]));
nexports_win = nexports(max_nexports<=prctile(v_netexports,cutoff),:);
nfa_t_win = nfa_t(max_nexports<=prctile(v_netexports,cutoff),:);

figure;histogram(nfa_t_win);title('histogram NFA');
figure;histogram(nexports_win);title('histogram Net exports');

errors_win = errors(max_nexports<=prctile(v_netexports,cutoff),:);
my_mean_errors = (mean(mean(errors_win)))
my_max_errors = max(max(errors_win))
