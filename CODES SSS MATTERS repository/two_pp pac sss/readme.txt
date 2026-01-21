These Matlab codes replicate the results from two-parameter perturbation around the stochastic steady state when calibrating the model with two-parameter perturbation (SSS solution), using Levintal's toolkit.

The 5 current scripts in this folder are:

fgru_pac_sss_errors: (called by solution_fgru_pac_sss_est) computes Euler equation errors from simulations (Table 1)

routine_sss_tasks: (called by various scripts) performs routine tasks 

solution_fgru_pac_sss_est: performs the estimation of parameters (Table 1)
solution_fgru_pac_sss_irf: computes IRFs (for figures 3, 4, 5, 6, 7, 9, 10, 13, and 14)
solution_fgru_pac_sss_var: computes simulated moments for variance decomposition (Table 2)

INSTRUCTIONS: 
- Add the parent folder of this repository to the Matlab path (folder only).
- Make sure you have run the script "prepare_model_fgru_pac_sss" (found in the folder "model pac sss") at least once  (make sure the folder is the main MATLAB path when executing)
