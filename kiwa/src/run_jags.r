#!/usr/bin/env Rscript --vanilla

# ==== Input parameters ====

'
Usage:
run_jags.r <dat> <mod> <ctf> <out> 
run_jags.r (-h | --help)

This script fits a Bayesian multi-state, robust design, Cormack-Jolly-Seber 
apparent survival model using the `jagsUI` interface to JAGS. User must supply a
pre-made .jags model file and prepared data.  Model results stored

Parameters:
dat: path to .Rdata file with data pre-prepared for JAGS
mod: model number
ctf: path to model control file
out: path to output directory

Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

# source("src/input_parse.r")

if(interactive()) {
  library(here)
  
  .wd <- '/home/sy522/project/kiwa_survival_analysis'
  
  rd <- here
  
  # TODO: check/fix filepaths
  .datPF <-file.path(.wd,'output/kiwa_dat_JAGS.rdata')
  .outPF <- file.path(.wd,'output')
  .ctfs <- file.path(.wd, 'ctfs/model_control.csv')
  
  mod <- 1
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  .script <-  thisfile()
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source(rd('src/input_parse.r'))
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datPF <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
  .ctfs <- makePath(ag$ctf)
  mod <- ag$mod
  
}

# ==== Setup ==== #

suppressWarnings(
  suppressPackageStartupMessages({
    library(jagsUI)
    library(tidyverse)
    library(glue)
  }))

#-- Load Data and CTFs --#
message("Loading data...")
load(.datPF)

message("Loading control files...")
ctf <- read.csv(.ctfs, stringsAsFactors = F) %>% 
  filter(model_no == mod)


# ==== Run Models ==== #

#- Declare model input

# bad hardcode for how many age classes for which we should generate inits.
age_classes <- 1

mean.phiW <- if(ctf$wPar==0){rep(0.99, age_classes)}else{NA}
mean.phiM <- if(ctf$mPar==0){rep(0.99, age_classes)}else{NA}
mean.phiB <- if(ctf$bPar==0){rep(0.99, age_classes)}else{NA}
muW <- if(ctf$wPar==0){NA}else{rep(0, age_classes)}
betaW <- if(ctf$wPar==0){NA}else{rep(0, age_classes)}
muM <- if(ctf$mPar==0){NA}else{rep(0, age_classes)}
betaM <- if(ctf$mPar==0){NA}else{rep(0, age_classes)}
muB <- if(ctf$bPar==0){NA}else{rep(0, age_classes)}
betaB <- if(ctf$bPar==0){NA}else{rep(0, age_classes)}

# Initial values
inits <- function(){
  l <- list(
    mean.phiW = mean.phiW,
    mean.phiM = mean.phiM,
    mean.phiB = mean.phiB,
    mean.pW = 0.25,
    mean.pM = 0.25, 
    mean.pB = 0.5,
    z = z.init,
    muW = muW,
    betaW = betaW,
    muM = muM,
    betaM = betaM,
    muB = muB,
    betaB = betaB
  )
  out <- l[!is.na(l)]
  return(out)
}

var <- if(is.na(ctf$var)){NA}else{get(ctf$var)}

## Bundle data
jags.data <- list(first = first,
                  ch = ch, 
                  n.ind = dim(ch)[1], 
                  n.prim = dim(ch)[2],
                  n.sec = n.secondary,
                  yesW = yesW,
                  yesM = yesM,
                  yesB = yesB,
                  totalW = totalW,
                  totalM = totalM,
                  totalB = totalB,
                  z = z,
                  x = ifelse(!is.na(var),var,0)
                  # age = age
                  # br = breeder

)


### Parameters
params <- c(
  "mean.phiW", "mean.phiM", "mean.phiB",
  "mean.pW", "mean.pM", "mean.pB",
  "gP_W", "gP_M", "gP_B",
  "gPP_W", "gPP_M", "gPP_B",
  "gWM", "kWM", "gMB", "kMB",
  "depWM", "depMB",
  "muW", "betaW",
  "muM", "betaM",
  "muB", "betaB"
  # "z"
  )

### MCMC settings
nc <- 4
ni <- 500000
nb <- 50000
nt <- 250
na <- 1000

# Fit model and view output
mod <- jagsUI::jags(data = jags.data, inits = inits, 
                    parameters.to.save = params, 
                    model.file = glue("jags/{ctf$jags_file}"), 
                    n.chains = nc, 
                    n.iter = ni, n.burnin = nb, n.thin = nt, 
                    n.adapt = na,
                    parallel = TRUE)

save(mod, file = glue("{.outPF}/model_w2017_noage_{ctf$model_no}.rdata"))
