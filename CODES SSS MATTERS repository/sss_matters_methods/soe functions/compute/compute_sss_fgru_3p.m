%This function approximates the SSS in the FGRU model, for the "SSS solution" case.
function y = compute_sss_fgru_3p(model,params,my_dss_params,M,eta,eps_ind,approx,guess)

my_evalf= @(x)eval_sss_fgru_3p(x,model,params,my_dss_params,M,eta,eps_ind,approx);
y = fsolve(my_evalf,guess); %,'Display','off'));