Brief description of each folder in this repository:

- estimated parameters: contains matlab .mat files
- sss_matters_methods: all functions and toolboxes used to replicate the main results (excluding Taylor projection)

- model pac dss: contains script to do pre-processing of the plain FGRU model, to be executed just once
- model pac loc: contains script to do pre-processing of the auxiliary perturbation model to approximate around the SSS ("SSS fixed" solutions), to be executed just once
- model pac sss: contains script to do pre-processing of the auxiliary perturbation model to approximate around the SSS (SSS solutions), to be executed just once
- model uza sss: contains script to do pre-processing of the auxiliary perturbation model (to compute both DSS and SSS solutions), to be executed just once


- one_pp pac dss: scripts to replicate results from standard perturbation around DSS
- one_pp uza dss: scripts to replicate results from standard perturbation around DSS (Uzawa model)
- two_pp pac loc: scripts to replicate results from two-parameter perturbation around SSS ("SSS fixed" solutions)
- two_pp pac sss: scripts to replicate results from two-parameter perturbation around SSS (SSS solutions)
- two_pp pac rec: scripts to replicate results from two-parameter perturbation around SSS (SSS "rec" solutions)
- two_pp uza sss: scripts to replicate results from two-parameter perturbation around SSS (Uzawa model)

- dynare codes: files to replicate impulse responses using Dynare (figures 2, 3 and 5 for DSS solution, figure 2 for SSS "fixed" solution)
- taylor projection codes: files and toolboxes to replicate results from Taylor projection

