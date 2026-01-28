Here's a description of the content in this folder:

- The folder "model pac dss" consists of files to solve the FGRU model around the DSS with the DSS calibration, and they replicate IRFs in figs. 2, 3 and 5.
These files are:
prepare_model_pac_dss.mod: Dynare's MOD file
pre_processing_pac_dss: MATLAB script, to execute the model's pre-processing
dyn_solution_fgru_pac_dss_irf: MATLAB script, computes the IRFs

- The folder "model pac loc" consists of files to solve the FGRU model around the SSS with the DSS calibration (SSS "fixed" solution), and they replicate IRFs in Fig 2.
These files are:
prepare_model_pac_loc.mod: Dynare's MOD file
pre_processing_pac_loc: MATLAB script, to execute the model's pre-processing
dyn_solution_fgru_pac_loc_irf: MATLAB script, computes the IRFs

In addition, the subfolder "functions" consists of the following Matlab functions:
compute_sss_fgru_3_locp_num_dyn: solves for the SSS, using the SSS algorithm described in Appendix C
eval_sss_fgru_3_locp_num_dyn: implements the SSS algorithm
my_dss_psi_k: auxiliary function to solve for the DSS
my_num_dss: computes DSS of auxiliary model numerically

INSTRUCTIONS: 
- Add the parent folder of this repository 'CODES SSS MATTERS repository' to the Matlab path (folder only).
- At the start of each Matlab section, execute the scripts 'pre_processing_pac_dss' and 'pre_processing_pac_loc' once, 
  making sure the Dynare installation (version 5.X) is in the correct path (by default, 'C:\dynare\5.5'), otherwise you must change this manually in codes.
NOTE: The Dynare toolbox and Levintal's don't get along too well.. So it is recommended to start a new MATLAB session if, after executing Dynare's codes, 
  you want to go back to work with the main files in the repository. Alternatively, you can try to remove Dynare's Matlab paths manually.