## ---------------------------
##
## Script name: buildjobs.r
##
## Purpose of script: Build bash script that will run run_jags.r duynamically to 
## fit a set of CJS models 
##
## Author(s): Shubhi Sharma, Scott Yanco (modified Shubi's original)
##
## Date Modified: 2021-10-01
##
## Email: scott.yanco@yale.edu
##
## ---------------------------


# wehn set to 1, this script will automatically submit the jobs generated
AUTO_SUBMIT <- 1


library(glue)

#filepaths

# This figures out if I'm on the HPC or local machine based on username 
user <- Sys.info()["user"]
if(user == "sy522"){
  project.directory <- "/home/sy522/project/kiwa_survival_analysis"
}else{
  project.directory <- "/home/sy522/project/kiwa_survival_analysis"
}

data.directory <- file.path(project.directory, "output")
scripts.directory <- file.path(project.directory, "src")
control.directory <- file.path(project.directory, "ctfs")
job.directory <- file.path(project.directory, "hpc_jobs")

ctf <- read.csv(file.path(control.directory, "model_control.csv"), 
                stringsAsFactors = F)

Rscript_prefix <- "Rscript " # generally "Rscript " will work, but in case you need to point to another install of Rscript, change this

mod_names <- ctf$formula
mod_nums <- ctf$model_no


header <- "module load miniconda
conda activate kiwa2"
shebang <- '#!/bin/bash'


prefix <- "#SBATCH --"

for(o in 1:length(mod_names)){
  jobname <- paste0(ctf$formula[o])
  job.file <- glue("job_model_{ctf$model_no[o]}.sh")
  
  command <- paste(Rscript_prefix,
                    file.path(scripts.directory, "run_jags.r"),
                    file.path(data.directory, "kiwa_dat_JAGS.rdata"),
                    ctf$model_no[o],
                    file.path(control.directory, "model_control.csv"),
                    data.directory,
                    sep = " ")
  
  mid <- paste0(prefix, "job-name=", jobname, "\n",
		prefix, "partition=day,pi_jetz\n",
                prefix, "time=2-00:00:00\n",
                prefix, "cpus-per-task=4\n",
                prefix, "mail-user=scott.yanco@yale.edu\n",
                prefix, "mail-type=ALL")
  
  # command_fin <- paste(command, models[o], dataset)
  
  writeLines(c(shebang, mid, "\n", header, command),
             file.path(job.directory, job.file))
  
  if(AUTO_SUBMIT & user == 'sy522')system(paste('sbatch',
                                                 file.path(job.directory, job.file)))
  if(!AUTO_SUBMIT){print(paste0("To submit, run
                                sbatch ", file.path(job.directory, job.file)))}
}
















