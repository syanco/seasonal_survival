# Analysis Workflow Associated with Cooper et al. *in prep*

Analyses for this paper are separated by the two focal species:  Kirtland's Warbler (KIWA) and Black-throated Blue Warbler (BTBW).

## Data Notes

Data associated with this project is store4d separately from this repository and can be found at: [PENDING]

Each analysis below describes the anticipated directory structure for input data.


## Kirtlands Warbler Analysis

This analysis relies on a series of scripts to execute the complete workflow.  Model fitting is expected to occur on an high-performance computing (HPC) environment (cluster) using a SLURM manager.  With modification, the code can be run without an HPC and notes below suggest how that can be accomplished.

### Data and Intermediate Products

Data should be stored in a directory named `/data` at the same level as, e.g., `src`. 

Similarly, the directory `/output` should be at the same level and is used to store derived outputs, including those needed as intermediary products.

Finally, the `/KIWA` directory should also include a folder named `/hpc_jobs` which stored .sh scripts generated programattical and used to launch model runs on the HPC environment.

For simplicity these three folders are included in this repository (but are empty) to facilitate cloning directly from github.

### Workflow

This section describes the "recipe" to reproduce the KIWA analysis.  Unless otherwise noted, scripts can be run either interactively or from the command line using `Rscript ...`

1.  Environmental Annotations:  `annotate_kiwa.r`  
*Annotates KIWA winter centroids with March EVI*

2.  Initial Data Processing: `prep_data_w2017.r`  
*Primary data processing script.  Links captures with environmental annoatations and individual morphology data.  Formats for survival analysis using JAGS and saves all objects needed to run the model as an .Rdata file*

3.  Write JAGS Models: `write_jags_noage.r`  
*Writes 8 JAGS scripts, one for each possible combination of seasons with covariates plus a model with no covariates. (Note that only the fully specified model and the "dot" model are considered in this analysis).  Files are saved as .jags format in the `/jags` directory.*

4.  Launch Models on HPC: `build_jobs.r`. 
*This script will automatically create .sh scripts that can be used to launch HPC jobs to run the JAGS models.  It will generate one script per covariate combination requested - these can be defined in* `/ctfs/model_control.csv`. *The .sh scripts will be stored in* `/hpc_jobs`. *If* `AUTO_SUBMIT` *variable is set to* `1` *at line 18 (the default), then the script will also send the* `sbatch` *command as welland will submit the job to the slurm.  Thus, we recommned running this from the HPC environment itself.  The .sh scripts ultimately call the* `run_jags.r` *script to actually draw MCMC samples.  If not using an HPC, this script can be used to manually run the JAGS model.  The* `run_jags.r` *script is set up to take command line options when called via* `Rscript ...` *or run interactively by modifying lines 35-39.*




## Black-throated Blue Warbler Analysis
