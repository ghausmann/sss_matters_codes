% This function implements the SSS algorithm as described in Appendix C.
%The key inputs are D1 (aux. debt level),lK1 (log of DSS capital), and
%psi_i (aux. depreciation parameter)

%The standard tasks are:
%1) using the key inputs, recalculate the DSS
%2) compute k-order derivatives of the decision rules
%3) evaluate the Taylor series imposing sigma=epsilon=1
%4) return h1(0,1,1).

%Other inputs are Dynares structs (M_, options_, and oo_), 'my_dss_params' (a
%subset of fixed parameter values), 'my_ind' (indices of variables in
%Dynare's vectors), and 'nshocks' (number of shocks).

function y = eval_sss_fgru_3_locp_num_dyn(x,M_, options_, oo_,my_dss_params,my_ind,nshocks)

eps_ind = my_ind(1); %index of perturbation param WITHIN state variables
st_ini_ind = my_ind(2); %index of first state variable
st_end_ind = my_ind(3); %index of last state variable
D_ind = my_ind(4); %index of debt state
K_ind = my_ind(5); %index of capital state
I_ind = my_ind(6); %index of investment state


D1 = x(1); % approximation point for debt
lK1 = x(2); % approximation point for capital
K1 = exp(lK1);
psi_i = x(3);
delta = M_.params(strcmp(M_.param_names,'delta')==1);

%Recompute DSS of auxiliary model
[yss,psi_k] = my_num_dss(x,my_dss_params);

M_.params(strcmp(M_.param_names,'Da')==1)    = D1;
M_.params(strcmp(M_.param_names,'psi_k')==1) = psi_k;
M_.params(strcmp(M_.param_names,'psi_i')==1) = psi_i;
oo_.steady_state = yss;

%Compute up to k-order derivatives
%use Dynare's function resol.m
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%STEP 3: Evaluate the k-order policy rule for states at the DSS,
%but imposing the model of interest (epsilon=sigma=1).
x0 = yss(st_ini_ind:st_end_ind,1); % start at the steady-state
x0(eps_ind) = 1; % impose the model of interest

%Evaluate the policy rule
%uses my own function dr_yt.m
y1 = dr_yt(mdr,yss,2,x0-yss(st_ini_ind:st_end_ind,1),zeros(nshocks,1));

%Return h1(0,1,1) 
eq1 = y1(D_ind) - D1;
eq2 = y1(K_ind) - (lK1);
eq3 = y1(I_ind) - log((delta+psi_i)*K1);

y = [eq1;eq2;eq3];