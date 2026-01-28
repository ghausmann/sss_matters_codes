These Matlab codes replicate the results from standard perturbation around the deterministic steady state (Uzawa model) using Levintal's toolkit

The 4 current scripts in this folder are:

fgru_uza_dss_errors: (called by solution_fgru_uza_dss_est) computes deterministic debt transition and Euler equation errors from simulations

routine_sol_dss_tasks_uza: (called by various scripts) performs routine tasks

solution_fgru_uza_dss_est: performs the estimation of parameters (Table 5)
solution_fgru_uza_dss_irf: computes IRFs (for Figure 12)

INSTRUCTIONS: 
- Add the parent folder of this repository 'CODES SSS MATTERS repository' to the Matlab path (folder only).
- Make sure you have run the script "prepare_model_fgru_uza_dss" (found in the folder "model uza sss") at least once (make sure the folder is the main MATLAB path when executing)
