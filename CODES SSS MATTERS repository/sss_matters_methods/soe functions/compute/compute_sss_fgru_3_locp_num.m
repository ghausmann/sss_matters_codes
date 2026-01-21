%This function approximates the SSS in the FGRU model, for the "SSS fixed" case.

function y = compute_sss_fgru_3_locp_num(model,params,my_dss_params,M,eta,eps_ind,approx,guess)
%This function approximates the SSS of the FGRU model in the "SSS fixed" case.

my_evalf= @(x)eval_sss_fgru_3_locp_num(x,model,params,my_dss_params,M,eta,eps_ind,approx);
%y = fsolve(my_evalf,guess); %,'Display','off'));
%options = optimoptions('fsolve','Tolfun',1e-9,'Display','off'); 
options = optimoptions('fsolve','Tolfun',1e-9); 
y = fsolve(my_evalf,guess,options); %,'Display','off'));