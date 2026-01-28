These Matlab codes replicate the results from two-parameter perturbation around the stochastic steady state when keeping the calibration from standard perturbation intact ("SSS fixed" solution), using Levintal's toolkit.

The 7 current scripts in this folder are:

fgru_pac_loc_errors: (called by solution_fgru_pac_loc_moments) computes Euler equation errors from simulations (Tables 1, 3, and 5)

routine_sss_loc_tasks: (called by various scripts) performs routine tasks 

solution_fgru_pac_loc_derivatives: computes consumption elasticities about the SSS (Table 4)
solution_fgru_pac_loc_derivatives_trans: replicates Figure 8
solution_fgru_pac_loc_irf: computes IRFs (for figures 2, 6, 7 and 11)
solution_fgru_pac_loc_moments: computes simulated moments (Tables 1, 3, and 5)
solution_fgru_pac_loc_var: computes simulated moments for variance decomposition (Table 2)

INSTRUCTIONS: 
- Add the parent folder of this repository 'CODES SSS MATTERS repository' to the Matlab path (folder only).
- Make sure you have run the script "prepare_model_fgru_pac_loc" (found in the folder "model pac loc") at least once (make sure the folder is the main MATLAB path when executing)
