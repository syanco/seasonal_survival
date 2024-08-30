# Analysis Workflow Associated with Cooper et al. 2024
Cooper, N.W., S.W. Yanco, Clark S. Rushing, T. Scott Sillett, and P.P. Marra. 2024 . Non-breeding conditions induce a carry-over effect on survival of migratory birds. Current Biology.

Analyses for this paper are separated by the two focal species:  Kirtland's Warbler (KIWA) and Black-throated Blue Warbler (BTBW).

The active and up to date version of this repository can be found at https://github.com/syanco/seasonal_survival/

## Data Notes

Data associated with this project is stored separately from this repository and can be found at: https://doi.org/10.25573/data.25436914

Each analysis below describes the anticipated directory structure for input data.


## Kirtland's Warbler Analysis

This analysis relies on a series of scripts to execute the complete KIWA workflow and can be found within the `/KIWA` directory.  Model fitting is expected to occur on a high-performance computing (HPC) environment (cluster) using a SLURM manager.  With modification, the code can be run without a HPC and notes below suggest how that can be accomplished.

### Data and Intermediate Products

Data should be stored in a directory named `/data` at the same level as, e.g., `src`. 

Similarly, the directory `/output` should be at the same level and is used to store derived outputs, including those needed as intermediary products.

Finally, the `/KIWA` directory should also include a folder named `/hpc_jobs` which stores the .sh scripts generated programattically. These are used to launch model runs on the HPC.

For simplicity these three folders are included in this repository (but are empty) to facilitate cloning directly from github.

### Workflow

This section describes the "recipe" to reproduce the KIWA analysis.  Unless otherwise noted, scripts can be run either interactively or from the command line using `Rscript ...`

1.  Environmental Annotations:  `annotate_kiwa.r`  
*Annotates KIWA winter centroids with March EVI*

2.  Initial Data Processing: `prep_data_w2017.r`  
*Primary data processing script.  Links captures with environmental annoatations and individual morphology data.  Formats for survival analysis using JAGS and saves all objects needed to run the model as an .Rdata file*

3.  Write JAGS Models: `write_jags_noage.r`  
*Writes 8 JAGS scripts, one for each possible combination of seasons with covariates plus a model with no covariates. (Note that only the fully specified model and the "dot" model are considered in this analysis).  Files are saved as .jags format in the `/jags` directory.*

4.  Launch Models on HPC: `build_jobs.r`  
*This script will automatically create .sh scripts that can be used to launch HPC jobs to run the JAGS models.  It will generate one script per covariate combination requested - these can be defined in* `/ctfs/model_control.csv`. *The .sh scripts will be stored in* `/hpc_jobs`. *If* `AUTO_SUBMIT` *variable is set to* `1` *at line 18 (the default), then the script will also send the* `sbatch` *command as welland will submit the job to the slurm.  Thus, we recommned running this from the HPC environment itself.  The .sh scripts ultimately call the* `run_jags.r` *script to actually draw MCMC samples.  If not using an HPC, this script can be used to manually run the JAGS model.  The* `run_jags.r` *script is set up to take command line options when called via* `Rscript ...` *or run interactively by modifying lines 35-39.*

5.  Summarize results: `results_scratch_noage.r`  
*This script evaluates information theoretic support amongst models and summarizes results from the top model*
**TODO:** This is still a scratch script that clould use some cleaning (and contians duplicative code for predicting conditional effects but it is at least complete).

6.  Make effects predictions: `make_conditional_predictions.r`  
*This script makes conditional effects predictions from the top (full) model.  These predictions are used to generate plots that depict the effects of the fitted model (in "Combined Processes", below).*

7. Mortality risk comparisons: `compare_risks.r`  
*Makes comparisons among seasons and EVI values to relativze risk both normalized per unit time and givenexpected duration of risk exposure*

8.  Population growth projections: `seasonal_matrix_model.r`  
*Creates a simple matrix projection model for populatioon growth and links winter, spring, and summer seasonal survival to EVI per the final fitted model in this analysis and is used to estimate lambda conditional on EVI.  Produces the plot used as Figure 4 in the paper.* 


## Black-throated Blue Warbler Analysis

The complete workflow for this species is contained in `BTBW_analysis.R` (except for odds ratio calculations below).

`btbw_mort_odds.R`  
*Makes comparisons among seasons and EVI values to relativize risk both normalized per unit time and givenexpected duration of risk exposure*



## Combined Processes
Steps in this section apply to both KIWA and BTBW datasets and primarily relate to plotting figures and summarizing results.


1.  Make conditional effects plots: `plot_conditionals.r`  
*Using conditional predictions generated by both the KIWA and BTBW workflows, this produces conditional effects plots that serve as the basis for FIgures 2 & 3 in the manuscript (less some markups added after the fact).*

