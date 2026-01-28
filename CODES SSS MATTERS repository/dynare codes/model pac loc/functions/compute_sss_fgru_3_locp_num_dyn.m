function y = compute_sss_fgru_3_locp_num_dyn(M_,options_,oo_,my_dss_params,my_ind,nshocks,guess)
% This function approximates the SSS of Home and Foreign bonds.

% Choose a second-order solution (enough for this model)
options_.order=2;
options_.k_order_solver=0;

%Call fsolve to find a candidate a_bar that zeroes the residual SSS
%function:
my_evalf= @(x)eval_sss_fgru_3_locp_num_dyn(x,M_, options_, oo_,my_dss_params,my_ind,nshocks);


%y = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));
y = fsolve(my_evalf,guess,optimset('MaxFunEval',10000));