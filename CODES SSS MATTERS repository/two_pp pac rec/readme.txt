These Matlab codes replicate the results from two-parameter perturbation around the stochastic steady state when recalibrating the model with two-parameter perturbation (SSS "rec" solution), using Levintal's toolkit.

The 5 current scripts in this folder are:

fgru_pac_rec_errors: (called by solution_fgru_pac_rec_est) computes Euler equation errors from simulations and erdogic distribution of external debt (Table 1 and Figure 1)

routine_rec_tasks: (called by various scripts) performs routine tasks 

solution_fgru_pac_rec_est: performs the estimation of parameters (Table 1)
solution_fgru_pac_rec_irf: computes IRFs (for figures 3, 4, and 5)
solution_fgru_pac_rec_var: computes simulated moments for variance decomposition (Table 2)

INSTRUCTIONS: 
- Add the parent folder of this repository 'CODES SSS MATTERS repository' to the Matlab path (folder only).
- Make sure you have run the script "prepare_model_fgru_pac_sss" (found in the folder "model pac sss") at least once  (make sure the folder is the main MATLAB path when executing)
