library(dplyr)
library(nimble)

############################
#### Nimble code ----
############################

n_code <- nimbleCode({
  ###################################################################################################
  ## 1. Priors
  ###################################################################################################
  
  ## ------------------------------------ 1.1 Seasonal survival ---------------------------------------
  ##-------------------------------
  ## 1.1a Season survival estimates
  ##-------------------------------
  lmean.phi.m.win <- logit(mean.phi.m.win)  # Logit mean male winter weekly survival probability
  mean.phi.m.win ~ dbeta(100, 1)       # Mean male weekly winter survival
  
  lmean.phi.m.sum <- logit(mean.phi.m.sum)  # Logit mean male summer weekly survival probability
  mean.phi.m.sum ~ dbeta(100, 1)       # Mean male weekly summer survival
  
  lmean.phi.f.win <- logit(mean.phi.f.win)  # Logit mean male winter weekly survival probability
  mean.phi.f.win ~ dbeta(100, 1)       # Mean male weekly winter survival
  
  lmean.phi.f.sum <- logit(mean.phi.f.sum)  # Logit mean male summer weekly survival probability
  mean.phi.f.sum ~ dbeta(100, 1)       # Mean male weekly summer survival
  
  ### Mean weekly survival probability rates for males
  lmean.phi.m.aut <- logit(mean.phi.m.aut)  # Logit mean male autumn weekly survival probability
  mean.phi.m.aut ~ dbeta(100, 1)         # Mean male autumn weekly survival probability
  
  lmean.phi.m.spr <- logit(mean.phi.m.spr)  # Logit mean male spring weekly survival probability
  mean.phi.m.spr ~ dbeta(100, 1)         # Mean male spring weekly survival probability
  
  
  ### Mean weekly survival probability rates for females
  lmean.phi.f.aut <- logit(mean.phi.f.aut)  # Logit mean female autumn weekly survival probability
  mean.phi.f.aut ~ dbeta(100, 1)         # Mean female autumn weekly survival probability
  
  lmean.phi.f.spr <- logit(mean.phi.f.spr)  # Logit mean female spring weekly survival probability
  mean.phi.f.spr ~ dbeta(100, 1)         # Mean female spring weekly survival probability
  
  
  ### Annual variation
  for (t in 1:nYears){
    z.phi.m.sum[t] ~ dnorm(0, sd = 1) # Annual variation in male summer weekly survival probability
    z.phi.f.sum[t] ~ dnorm(0, sd = 1) # Annual variation in female summer weekly survival probability
    
    z.phi.m.win[t] ~ dnorm(0, sd = 1) # Annual variation in male winter weekly survival probability
    z.phi.f.win[t] ~ dnorm(0, sd = 1) # Annual variation in female winter weekly survival probability
    
    z.phi.m.spr[t] ~ dnorm(0, sd = 1) # Annual variation in male spring weekly survival probability
    z.phi.f.spr[t] ~ dnorm(0, sd = 1) # Annual variation in female spring weekly survival probability
  }
  
  for (t in 1:(nYears - 1)){
    z.phi.m.aut[t] ~ dnorm(0, sd = 1) # Annual variation in male autumn weekly survival probability
    z.phi.f.aut[t] ~ dnorm(0, sd = 1) # Annual variation in female autumn weekly survival probability
  }
  
  sigma.phi.sum ~ dgamma(1, 1) # Standard deviation of annual variation in summer weekly survival probability
  
  sigma.phi.win ~ dgamma(1, 1) # Standard deviation of annual variation in winter weekly survival probability
  
  sigma.phi.spr ~ dgamma(1, 1) # Standard deviation of annual variation in spring weekly survival probability
  
  sigma.phi.aut ~ dgamma(1, 1) # Standard deviation of annual variation in autumn weekly survival probability
  
  ### Slope coefficients
  beta.m.soi[1] ~ dnorm(0, sd = 1.75) # Male, winter
  beta.m.soi[2] ~ dnorm(0, sd = 1.75) # Male, spring
  beta.m.soi[3] ~ dnorm(0, sd = 1.75) # Male, breeding
  
  beta.f.soi[1] ~ dnorm(0, sd = 1.75) # Female, winter
  beta.f.soi[2] ~ dnorm(0, sd = 1.75) # Female, spring
  beta.f.soi[3] ~ dnorm(0, sd = 1.75) # Female, breeding
  
  beta.m.bseas[1] ~ dnorm(0, sd = 1.75) # Male, breeding
  beta.m.bseas[2] ~ dnorm(0, sd = 1.75) # Male, autumn
  beta.m.bseas[3] ~ dnorm(0, sd = 1.75) # Male, winter
  
  beta.f.bseas[1] ~ dnorm(0, sd = 1.75) # Female, breeding
  beta.f.bseas[2] ~ dnorm(0, sd = 1.75) # Female, autumn
  beta.f.bseas[3] ~ dnorm(0, sd = 1.75) # Female, winter
  
  ##-----------------------------
  ## 1.2b Detection probabilities
  ##-----------------------------
  
  ### Mean
  lmean.p.may <- logit(mean.p.may)  # Logit mean male May detection probability
  mean.p.may ~ dbeta(1, 1)            # Mean male May detection probability
  
  lmean.p.m.aug <- logit(mean.p.m.aug)   # Logit mean male August detection probability
  mean.p.m.aug ~ dbeta(1, 1)            # Mean male August detection probability 
  
  lmean.p.f.aug <- logit(mean.p.f.aug)   # Logit mean female August detection probability
  mean.p.f.aug ~ dbeta(1, 1)            # Mean female August detection probability 
  
  lmean.p.oct <- logit(mean.p.oct)     # Logit mean October detection probability
  mean.p.oct ~ dbeta(1, 1)             # Mean October detection probability 
  
  lmean.p.m.mar <- logit(mean.p.m.mar)     # Logit mean March detection probability
  mean.p.m.mar ~ dbeta(1, 1)             # Mean  March detection probability 

  lmean.p.f.mar <- logit(mean.p.f.mar)     # Logit mean March detection probability
  mean.p.f.mar ~ dbeta(1, 1)             # Mean  March detection probability 
  
  
  ### Annual variation in detection probability
  for(t in 1:nYears){
    # Summer
    z.p.aug[t] ~ dnorm(0, sd = 1)

    # Winter
    z.p.mar[t] ~ dnorm(0, sd = 1)
  }
  
  sigma.p.sum ~  dgamma(1, 1)
  sigma.p.win ~  dgamma(1, 1)
  
  ###################################################################################################
  ## 2. Constrain parameters
  ###################################################################################################
  
  
  ## -------------------------------------- 2.1 Seasonal Survival -----------------------------------
  
  for (t in 1:nYears){
    # Over-winter survival
    logit(wphi.m.win[t]) <- lmean.phi.m.win + beta.m.soi[1] * soi[t] + beta.m.bseas[3] * bseas[t] + z.phi.m.win[t] * sigma.phi.win          
    logit(wphi.f.win[t]) <- lmean.phi.f.win + beta.f.soi[1] * soi[t] + beta.f.bseas[3] * bseas[t] + z.phi.f.win[t] * sigma.phi.win
    PHI.m.win[t] <- pow(wphi.m.win[t], nWeeks.win) 
    PHI.f.win[t] <- pow(wphi.f.win[t], nWeeks.win) 
    
    # Spring survival
    logit(wphi.m.spr[t])   <- lmean.phi.m.spr + beta.m.soi[2] * soi[t] + z.phi.m.spr[t] * sigma.phi.spr
    logit(wphi.f.spr[t])   <- lmean.phi.f.spr + beta.f.soi[2] * soi[t] + z.phi.f.spr[t] * sigma.phi.spr
    PHI.m.spr[t] <- pow(wphi.m.spr[t], nWeeks.spr) 
    PHI.f.spr[t] <- pow(wphi.f.spr[t], nWeeks.spr) 
    
    # Summer survival
    logit(wphi.m.sum[t]) <- lmean.phi.m.sum + beta.m.soi[3] * soi[t] + beta.m.bseas[1] * bseas[t + 1] + z.phi.m.sum[t] * sigma.phi.sum
    logit(wphi.f.sum[t]) <- lmean.phi.f.sum + beta.f.soi[3] * soi[t] + beta.f.bseas[1] * bseas[t + 1] + z.phi.f.sum[t] * sigma.phi.sum
    PHI.m.sum[t] <- pow(wphi.m.sum[t], nWeeks.sum) 
    PHI.f.sum[t] <- pow(wphi.f.sum[t], nWeeks.sum) 
  } # t
  
  
  # Autumn survival
  for (t in 1:(nYears - 1)){
    logit(wphi.m.aut[t]) <- lmean.phi.m.aut + beta.m.bseas[2] * bseas[t + 1] + z.phi.m.aut[t] * sigma.phi.aut
    logit(wphi.f.aut[t]) <- lmean.phi.f.aut + beta.f.bseas[2] * bseas[t + 1] + z.phi.f.aut[t] * sigma.phi.aut
    PHI.m.aut[t] <- pow(wphi.m.aut[t], nWeeks.aut) 
    PHI.f.aut[t] <- pow(wphi.f.aut[t], nWeeks.aut) 
  } # t
  
  
  
  ## --------------------------- 2.2 Seasonal detection probabilities ------------------------------
  for (t in 1:(nYears - 1)){ 
    logit(p.m.sum[2 * t - 1]) <- lmean.p.m.aug + z.p.aug[t] * sigma.p.sum
    logit(p.m.win[2 * t - 1]) <- lmean.p.m.mar + z.p.mar[t] * sigma.p.win
    logit(p.m.sum[2 * t]) <- lmean.p.may
    logit(p.m.win[2 * t]) <- lmean.p.oct
    
    logit(p.f.sum[2 * t - 1]) <- lmean.p.f.aug + z.p.aug[t] * sigma.p.sum
    logit(p.f.win[2 * t - 1]) <- lmean.p.f.mar + z.p.mar[t] * sigma.p.win
    logit(p.f.sum[2 * t]) <- lmean.p.may
    logit(p.f.win[2 * t]) <- lmean.p.oct
  } # t
  
  logit(p.m.sum[nYears * nSeasons - 1]) <- lmean.p.m.aug  +  z.p.aug[nYears] * sigma.p.sum
  logit(p.m.win[nYears * nSeasons - 1]) <- lmean.p.m.mar  +  z.p.mar[nYears] * sigma.p.win
  logit(p.f.sum[nYears * nSeasons - 1]) <- lmean.p.f.aug  +  z.p.aug[nYears] * sigma.p.sum
  logit(p.f.win[nYears * nSeasons - 1]) <- lmean.p.f.mar  +  z.p.mar[nYears] * sigma.p.win
  
  
  ##################################################################################################
  ## 3. Vectorize survival estimates for CJS models
  ##################################################################################################
  
  ## ------------------------- 3.1 Breeding survival estimates from New Hampshire ------------------
  for (t in 1:(nYears - 1)){
    phi.m.sum[2 * t - 1] <- PHI.m.sum[t]                         # Male breeding season survival
    phi.f.sum[2 * t - 1] <- PHI.f.sum[t]                         # Female breeding season survival
    
    phi.m.sum[2 * t]  <- PHI.m.aut[t] * PHI.m.win[t + 1] * PHI.m.spr[t + 1]  # Male between-breeding survival
    phi.f.sum[2 * t]  <- PHI.f.aut[t] * PHI.f.win[t + 1] * PHI.f.spr[t + 1]  # Female between-breeding survival
  } # t
  
  phi.m.sum[2 * nYears - 1] <- PHI.m.sum[nYears]
  phi.f.sum[2 * nYears - 1] <- PHI.f.sum[nYears]
  
  
  
  ## --------------------------- 3.2 Winter survival estimates from Jamaica ------------------------
  for (t in 1:(nYears - 1)){
    phi.m.win[2 * t - 1] <- PHI.m.win[t]
    phi.f.win[2 * t - 1] <- PHI.f.win[t]
    
    phi.m.win[2 * t]  <- PHI.m.spr[t] * PHI.m.sum[t] * PHI.m.aut[t]  # Male between-breeding survival
    phi.f.win[2 * t]  <- PHI.f.spr[t] * PHI.f.sum[t] * PHI.f.aut[t]  # Female between-breeding survival
  } # t
  
  phi.m.win[2 * nYears - 1] <- PHI.m.win[nYears]
  phi.f.win[2 * nYears - 1] <- PHI.f.win[nYears]
  
  ##################################################################################################
  ## 4. Likelihoods
  ##################################################################################################
  
  
  ## 4.1 Likelihood for summer capture-recapture data: m-array CJS model
  # Multinomial likelihood
  for(t in 1:(nSeasons * nYears - 1)){
    q.m.sum[t] <- 1 - p.m.sum[t]              # Probability of non-recapture for males
    q.f.sum[t] <- 1 - p.f.sum[t]              # Probability of non-recapture for females
  }
  
  for (t in 1:(nSeasons * nYears - 1)){
    marr.m.sum[t, 1:(nSeasons * nYears)] ~ dmulti(pr.m.sum[t, 1:(nSeasons * nYears)], r.m.sum[t])
    marr.f.sum[t, 1:(nSeasons * nYears)] ~ dmulti(pr.f.sum[t, 1:(nSeasons * nYears)], r.f.sum[t])
  }
  
  # Define the cell probabilities of the m-arrays
  # Main diagonal
  for (t in 1:(nSeasons * nYears - 1)){
    pr.m.sum[t, t]  <- phi.m.sum[t] * p.m.sum[t]
    pr.f.sum[t, t]  <- phi.f.sum[t] * p.f.sum[t]
    
    # Above main diagonal   
    for (j in (t + 1):(nSeasons * nYears - 1)){
      pr.m.sum[t, j] <- prod(phi.m.sum[t:j]) * prod(q.m.sum[t:(j - 1)]) * p.m.sum[j]
      pr.f.sum[t, j] <- prod(phi.f.sum[t:j]) * prod(q.f.sum[t:(j - 1)]) * p.f.sum[j]
    } #j
    
    # Below main diagonal
    for (j in 1:(t - 1)){
      pr.m.sum[t, j] <- 0
      pr.f.sum[t, j] <- 0
    } #j
  } #t
  
  # Last column: probability of non-recapture
  for (t in 1:(nSeasons * nYears - 1)){
    pr.m.sum[t, nSeasons * nYears] <- 1 - sum(pr.m.sum[t, 1:(nSeasons * nYears - 1)])
    pr.f.sum[t, nSeasons * nYears] <- 1 - sum(pr.f.sum[t, 1:(nSeasons * nYears - 1)])
  } #t
  
  ## 4.1b Posterior-predictive checks for summer capture-recapture data
  # Compute fit statistics for observed data
  for (t in 1:(nSeasons * nYears - 1)){
    for (j in 1:(nSeasons * nYears)){
      expmarr.m.sum[t, j] <- r.m.sum[t] * pr.m.sum[t, j]
      E.org.m.sum[t, j] <- pow((pow(marr.m.sum[t, j], 0.5) - pow(expmarr.m.sum[t, j], 0.5)), 2)
      
      expmarr.f.sum[t, j] <- r.f.sum[t] * pr.f.sum[t, j]
      E.org.f.sum[t, j] <- pow((pow(marr.f.sum[t, j], 0.5) - pow(expmarr.f.sum[t, j], 0.5)), 2)
    } #j
  } #t
  
  # Generate replicate data and compute fit stats from them
  for (t in 1:(nSeasons * nYears - 1)){
    marr.m.sum.new[t, 1:(nSeasons * nYears)] ~ dmulti(pr.m.sum[t, 1:(nSeasons * nYears)], r.m.sum[t])
    marr.f.sum.new[t, 1:(nSeasons * nYears)] ~ dmulti(pr.f.sum[t, 1:(nSeasons * nYears)], r.f.sum[t])
    for (j in 1:(nSeasons * nYears)){
      E.new.m.sum[t,j] <- pow((pow(marr.m.sum.new[t, j], 0.5) - pow(expmarr.m.sum[t, j], 0.5)), 2)
      E.new.f.sum[t,j] <- pow((pow(marr.f.sum.new[t, j], 0.5) - pow(expmarr.f.sum[t, j], 0.5)), 2)
    } #j
  } #t
  
  fit.m.sum <- sum(E.org.m.sum[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  fit.m.sum.new <- sum(E.new.m.sum[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  
  fit.f.sum <- sum(E.org.f.sum[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  fit.f.sum.new <- sum(E.new.f.sum[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  
  ## 4.2a Likelihood for winter capture-recapture data: m-array CJS model
  # Multinomial likelihood
  for(t in 1:(nSeasons * nYears - 1)){
    q.m.win[t] <- 1 - p.m.win[t]              # Probability of non-recapture for males
    q.f.win[t] <- 1 - p.f.win[t]              # Probability of non-recapture for female
  }
  
  for (t in 1:(nSeasons * nYears - 1)){
    marr.m.win[t, 1:(nSeasons * nYears)] ~ dmulti(pr.m.win[t, 1:(nSeasons * nYears)], r.m.win[t])
    marr.f.win[t, 1:(nSeasons * nYears)] ~ dmulti(pr.f.win[t, 1:(nSeasons * nYears)], r.f.win[t])
  }
  
  # Define the cell probabilities of the m-arrays
  # Main diagonal
  for (t in 1:(nSeasons * nYears - 1)){
    pr.m.win[t, t]  <- phi.m.win[t] * p.m.win[t]
    pr.f.win[t, t]  <- phi.f.win[t] * p.f.win[t]
    
    # Above main diagonal   
    for (j in (t + 1):(nSeasons * nYears - 1)){
      pr.m.win[t, j] <- prod(phi.m.win[t:j]) * prod(q.m.win[t:(j - 1)]) * p.m.win[j]
      pr.f.win[t, j] <- prod(phi.f.win[t:j]) * prod(q.f.win[t:(j - 1)]) * p.f.win[j]
      
    } #j
    
    # Below main diagonal
    for (j in 1:(t - 1)){
      pr.m.win[t, j] <- 0
      pr.f.win[t, j] <- 0
    } #j
  } #t
  
  # Last column: probability of non-recapture
  for (t in 1:(nSeasons * nYears - 1)){
    pr.m.win[t, nSeasons * nYears] <- 1 - sum(pr.m.win[t, 1:(nSeasons * nYears - 1)])
    pr.f.win[t, nSeasons * nYears] <- 1 - sum(pr.f.win[t, 1:(nSeasons * nYears - 1)])
  } #t
  
  ## 4.2b Posterior-predictive checks for wintercapture-recapture data
  # Compute fit statistics for observed data
  for (t in 1:(nSeasons * nYears - 1)){
    for (j in 1:(nSeasons * nYears)){
      expmarr.m.win[t, j] <- r.m.win[t] * pr.m.win[t, j]
      E.org.m.win[t, j] <- pow((pow(marr.m.win[t, j], 0.5) - pow(expmarr.m.win[t, j], 0.5)), 2)
      
      expmarr.f.win[t, j] <- r.f.win[t] * pr.f.win[t, j]
      E.org.f.win[t, j] <- pow((pow(marr.f.win[t, j], 0.5) - pow(expmarr.f.win[t, j], 0.5)), 2)
    } #j
  } #t
  
  # Generate replicate data and compute fit stats from them
  for (t in 1:(nSeasons * nYears - 1)){
    marr.m.win.new[t, 1:(nSeasons * nYears)] ~ dmulti(pr.m.win[t, 1:(nSeasons * nYears)], r.m.win[t])
    marr.f.win.new[t, 1:(nSeasons * nYears)] ~ dmulti(pr.f.win[t, 1:(nSeasons * nYears)], r.f.win[t])
    for (j in 1:(nSeasons * nYears)){
      E.new.m.win[t,j] <- pow((pow(marr.m.win.new[t, j], 0.5) - pow(expmarr.m.win[t, j], 0.5)), 2)
      E.new.f.win[t,j] <- pow((pow(marr.f.win.new[t, j], 0.5) - pow(expmarr.f.win[t, j], 0.5)), 2)
    } #j
  } #t
  
  fit.m.win <- sum(E.org.m.win[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  fit.m.win.new <- sum(E.new.m.win[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  
  fit.f.win <- sum(E.org.f.win[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  fit.f.win.new <- sum(E.new.f.win[1:(nSeasons * nYears - 1), 1:(nSeasons * nYears)])
  
  
  ##################################################################################################
  ## 5. Derived parameters
  ##################################################################################################
  
  for(t in 1:(nYears - 1)){
    PHI.m.ann[t] <- PHI.m.sum[t] * PHI.m.aut[t] * PHI.m.win[t] * PHI.m.aut[t]
    PHI.f.ann[t] <- PHI.f.sum[t] * PHI.f.aut[t] * PHI.f.win[t] * PHI.f.aut[t]
  }
})

############################
#### Nimble model ----
############################

### Load data (assumes btbw_data.Rdata has been downloaded. Change path to match file location)
load("data/btbw_data.RData")

n_model <- nimbleModel(code = n_code,
                       constants = list(nYears = nYears, 
                                        nSeasons = 2,
                                        nWeeks.aut = 8,
                                        nWeeks.spr = 8,
                                        nWeeks.win = 24,
                                        nWeeks.sum = 12),
                       inits = list(mean.phi.m.sum = runif(1, 0.99, 0.999), 
                                    mean.phi.m.win = runif(1, 0.99, 0.999),
                                    mean.phi.f.sum = runif(1, 0.99, 0.999), 
                                    mean.phi.f.win = runif(1, 0.99, 0.999),
                                    ## Male weekly survival probability
                                    mean.phi.m.aut = runif(1, 0.99, 0.999), 
                                    mean.phi.m.spr = runif(1, 0.99, 0.999), 
                                    ## Female weekly survival probability
                                    mean.phi.f.aut = runif(1, 0.99, 0.999),
                                    mean.phi.f.spr = runif(1, 0.99, 0.999),
                                    ## Annual variation in weekly survival probability
                                    z.phi.m.sum = rnorm(nYears), 
                                    z.phi.m.win = rnorm(nYears), 
                                    z.phi.m.aut = rnorm(nYears - 1),
                                    z.phi.m.spr = rnorm(nYears), 
                                    z.phi.f.sum = rnorm(nYears),
                                    z.phi.f.win = rnorm(nYears),
                                    z.phi.f.aut = rnorm(nYears - 1),
                                    z.phi.f.spr = rnorm(nYears),
                                    ## SD weekly survival probability
                                    sigma.phi.sum = runif(1), 
                                    sigma.phi.win = runif(1),
                                    sigma.phi.aut = runif(1), 
                                    sigma.phi.spr = runif(1), 
                                    ## Mean p
                                    mean.p.may = runif(1), 
                                    mean.p.m.aug = runif(1),
                                    mean.p.f.aug = runif(1),
                                    mean.p.oct = runif(1), 
                                    mean.p.m.mar = runif(1),
                                    mean.p.f.mar = runif(1),
                                    ## Annual variation in p
                                    z.p.mar = rnorm(nYears),
                                    z.p.aug = rnorm(nYears),
                                    ## SD p
                                    sigma.p.sum = runif(1), 
                                    sigma.p.win = runif(1),
                                    ## Slope coefficients
                                    beta.m.soi = rnorm(3, 0, 0.25), beta.m.bseas = rnorm(3, 0, 0.25),
                                    beta.f.soi = rnorm(3, 0, 0.25), beta.f.bseas = rnorm(3, 0, 0.25)), 
                       data = list(marr.m.sum = nh.m.marray, r.m.sum = r.nh.m,
                                   marr.f.sum = nh.f.marray, r.f.sum = r.nh.f,
                                   marr.m.win = jam.m.marray, r.m.win = r.jam.m,  
                                   marr.f.win = jam.f.marray, r.f.win = r.jam.f,
                                   bseas = stn.bseas, soi = soi))

n_model$initializeInfo()

##################################################
#### Configure model, build MCMC, compile ----
##################################################

n_config <- configureMCMC(n_model)

n_config$removeSamplers(c('mean.phi.m.aut', 'mean.phi.m.spr', 'mean.phi.f.aut', 'mean.phi.f.spr'))
n_config$addSampler(target = c('mean.phi.f.aut', 'mean.phi.f.spr'), type = 'AF_slice')
n_config$addSampler(target = c('mean.phi.m.aut', 'mean.phi.m.spr'), type = 'AF_slice')

n_config$addMonitors(c("PHI.m.sum", "PHI.m.win", "PHI.m.spr", "PHI.m.aut",
                       "PHI.f.sum", "PHI.f.win", "PHI.f.spr", "PHI.f.aut",
                       "PHI.m.ann", "PHI.f.ann",
                       "fit.m.sum", "fit.m.sum.new", "fit.f.sum", "fit.f.sum.new",
                       "fit.m.win", "fit.m.win.new", "fit.f.win", "fit.f.win.new"))

n_mcmc <- buildMCMC(n_config)

comp_mod <- compileNimble(n_model)
comp_mcmc <- compileNimble(n_mcmc)

#####################################
#### Fit model ----
#####################################

system.time(n_samples <- runMCMC(comp_mcmc,
                                 niter = 50000,
                                 nburnin = 5000,
                                 thin = 1, 
                                 nchains = 3,
                                 samplesAsCodaMCMC = TRUE))


MCMCvis::MCMCdiag(n_samples, 
         mkdir = paste0("Prob-", Sys.Date(), "-SillettHolmes"),
         file_name = paste0("Prob-", Sys.Date()),
         dir = 'results', 
         save_obj = TRUE,
         obj_name = 'model-output',
         add_obj = list(n_code, sessionInfo()),
         add_obj_names = c('model-code', 
                           'session-info'),
         cp_file = c('scripts/02-Analysis_probs.R'))
