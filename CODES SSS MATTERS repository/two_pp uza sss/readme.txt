These Matlab codes replicate the results from two-parameter perturbation around the stochastic steady state when calibrating the model with two-parameter perturbation (Uzawa model), using Levintal's toolkit.

The 4 current scripts in this folder are:

fgru_uza_sss_errors: (called by solution_fgru_uza_sss_est) computes Euler equation errors from simulations

routine_sss_tasks: (called by other scripts) performs routine tasks 

solution_fgru_pac_sss_est: performs the estimation of parameters (Table 5)
solution_fgru_pac_sss_irf: computes IRFs (for Figure 12)


INSTRUCTIONS: 
- Add the parent folder of this repository 'CODES SSS MATTERS repository' to the Matlab path (folder only).
- Make sure you have run the script "prepare_model_fgru_uza_sss" (found in the folder "model uza sss") at least once (make sure the folder is the main MATLAB path when executing)
