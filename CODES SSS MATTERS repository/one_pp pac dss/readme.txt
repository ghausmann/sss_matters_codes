These Matlab codes replicate the results from standard perturbation around the deterministic steady state using Levintal's toolkit

The 13 current scripts in this folder are:

fgru_pac_dss_errors: (called by solution_fgru_pac_dss_est and solution_fgru_pac_dss_est_checks) computes deterministic debt transition (Figure 1, panel (a)),
	             Euler equation errors from simulations, and ergodic distribution of external debt (Figure 1, panel (c))
fgru_pac_dss_errors_pctile: (called by solution_fgru_pac_dss_est_pctile) computes Euler errors from simulations for higher-order solutions (Table 3)

routine_est_dss_tasks: (called by various scripts) performs routine tasks 
routine_sol_dss_tasks: (called by various scripts) performs routine tasks

solution_fgru_pac_dss_derivatives: computes consumption elasticities about the DSS (Table 4)
solution_fgru_pac_dss_est: performs the estimation of parameters (Table 1)
solution_fgru_pac_dss_est_checks: performs the estimation of parameters for robustness checks (Table 5)
solution_fgru_pac_dss_est_pctile: computes simulated moments for higher-order solutions (Table 3)

solution_fgru_pac_dss_exports: replicates Figure 1 panels (b)

solution_fgru_pac_dss_irf: computes IRFs (for figures 2, 3, 5, 6, 7, 9 and 10)
solution_fgru_pac_dss_irf_checks: computes IRFs for robustness checks (for Figure 11)

solution_fgru_pac_dss_simul: simulation of exogenous states, so we can feed Taylor projection's simulations with them.

solution_fgru_pac_dss_var: computes simulated moments for variance decomposition (Table 2)

INSTRUCTIONS: 
- Add the parent folder of this repository to the Matlab path (folder only).
- Make sure you have run the script "prepare_model_fgru_pac_dss" (found in the folder "model pac dss") at least once (make sure the folder is the main MATLAB path when executing)
