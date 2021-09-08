# MRC-BSU-Internship
The aim of this work is to develop an extension to the Partial Ordering Continual Reassessment Method (POCRM) proposed by Wages et al. (2011) to take into account further uncertainty in the ordering of doses in settings where this uncertainty is present. This work was completed at the MRC Biostatistics Unit between June 2021 and September 2021. The documents and R code associated with this project, that can be found within this repository, are outlined below.

# Documents
Relevant presentations and documents produced throughout the internship. 
1. **BCRM.pdf:** *Presentation on the Bayesian Continual Reassessment Method*
2. **BMA_POCRM_Spec.pdf:** *Specification of the final implementations of the POCRM found in this work.*
3. **BMA_prior_elic.pdf:** *Outlines the basic point estimate method for BMA POCRM and introduces a novel prior elicitation method.*
4. **Dose_schedule_finding_trials.pdf:** *Outlines an extension of Jaki et al. (2021) to handle partial orderings in dose schedule finding trial.*
5. **equivalence_proof.pdf:** *Proves the equivalence of the basic point estimate and mixture distribution point estimate methods.*
6. **Geo_mean_pcs_delta.png:** *Graph of geometric mean of proportion of correct selections for various delta under each method.*
7. **geom_mean_pcs_delta_large.png:** *Graph of geometric mean of proportion of correct selections for various delta under each method.*
8. **NIHR_Methodology_Internship_Final_Presentation.pdf:** *Final presentation outlining all work completed.*

# R 
1. **bma_pocrm_post.R:** *Implementation of all methods outlined in BMA_POCRM_Spec.pdf.*
2. **bma_pocrm_test.R:** *Debugging implementations of all methods for single trial dose allocation functions.*
3. **bma_sim_test.R:** *Debugging implementations of all methods for repeated trial dose allocation functions.*
4. **bayesian_pocrm_sim.R:** *Replication of simulation studies found in Wages et al. (2011).*
5. **POCRM_sim_rep.Rmd:** *Results of replication of simulation studies found in Wages et al. (2011).*
6. **bma_single_sim.Rmd:** *Assessment of behaviour of all methods under various scenarios for single trial runs.*
7. **pocrm_mod.Rmd:** *Outline of modifications made to original POCRM.*
8. **pocrm_calib.R:** *Running parameter calibration via High Performance Computing (HPC) resources .*
9. **pocrm_sim.R:** *Running simulation study based on Barnett et al. (2021) via HPC resources.*
10. **calibration_analysis.R:** *Analyse parameter calibration results.*
11. **simulation_analysis.R:** *Analyse simulation study results.*
12. **slurm_calib.txt:** *Script used for calibration via HPC.*
13. **slurm_sim.txt:** *Script used for simulation study via HPC.*
14. **output:** *Contains output from calibration and simulation studies. Note: For analysis to be carried out, these files must be in the same folder as the corresponding analysis scripts.*
