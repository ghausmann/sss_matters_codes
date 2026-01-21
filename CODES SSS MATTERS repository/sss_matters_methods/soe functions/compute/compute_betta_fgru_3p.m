%This function approximates the SSS in the FGRU model, for the "SSS
%solution" case in "calibration" mode.
function y = compute_betta_fgru_3p(D1,model,params,my_dss_params,M,eta,eps_ind,approx,guess)

my_evalf= @(x)calib_sss_fgru_3p(x,D1,model,params,my_dss_params,M,eta,eps_ind,approx);
%keyboard;
y = fsolve(my_evalf,guess,optimset('TolFun',1e-10,'Display','off'));
%y = fsolve(my_evalf,guess,optimset('TolFun',1e-10));