####################################
####                            ####
####      Write Jags Models     ####
####      Scott Yanco, PhD      ####
####    scott.yanco@yale.edu    ####
####                            ####
####################################

##---- Description ----##

# This script writes the base jags models to file for Kirtland's Warbler (KIWA) 
# survival analysis.  Basic model is a multi-state Cormack-Jolly-Seber (CJS) where
# states are seasons (winter, migration, breeding). Model 
# variations can account for the inclusion of up to 1 covariate in each season. 


##---- Initialization ----##

# NA

##---- Write Models ----##

####-- "Dot" Model --####

sink("jags/ms-cjs_dot.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  
        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions
            # Winter and Migration vary by age only, Breeding varies by age and
            # breeder status.  No estimates vary by occasion (loop included for
            # future flexibility)
            
            phiW[  i, t] <- mean.phiW
            phiM[  i, t] <- mean.phiM
            phiB[  i, t] <- mean.phiB
          } #t
        } #i

    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #   depMB[i] <- 1
    # }
    
    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    
    # TODO: currently rpoducing impossible beta distributions, fix
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)
    

      mean.phiM ~ dbeta(3, 0.5)T(0.0001,0.9999)
      mean.phiW ~ dbeta(5, 0.5)T(0.0001,0.9999)
      mean.phiB ~ dbeta(5, 0.5)T(0.0001,0.9999)

    
    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)
    
    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)

    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)

  
  
  #-- Transition Matrices --#
  
      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1,i] <- phiW[i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2,i] <- phiW[i,t] * gPP_M * depWM[t]
          ps[1,t,3,i] <- 0
          ps[1,t,4,i] <- phiW[i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5,i] <- phiW[i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6,i] <- 0
          ps[1,t,7,i] <- (1 - phiW[i,t])
          
          ps[2,t,1,i] <- 0
          ps[2,t,2,i] <- phiM[i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3,i] <- phiM[i,t] * gPP_B * depMB[t] 
          ps[2,t,4,i] <- 0
          ps[2,t,5,i] <- phiM[i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6,i] <- phiM[i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7,i] <- (1 - phiM[i,t])
          
          ps[3,t,1,i] <- 0
          ps[3,t,2,i] <- 0 
          ps[3,t,3,i] <- phiB[i,t] * gPP_B
          ps[3,t,4,i] <- 0 
          ps[3,t,5,i] <- 0
          ps[3,t,6,i] <- phiB[i,t] * (1 - gPP_B)
          ps[3,t,7,i] <- (1 - phiB[i,t])
          
          ps[4,t,1,i] <- phiW[i,t] * gP_W * (1-depWM[t])
          ps[4,t,2,i] <- phiW[i,t] * gP_M  * depWM[t]
          ps[4,t,3,i] <- 0
          ps[4,t,4,i] <- phiW[i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5,i] <- phiW[i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6,i] <- 0
          ps[4,t,7,i] <- (1 - phiW[i,t]) 
          
          ps[5,t,1,i] <- 0
          ps[5,t,2,i] <- phiM[i,t] * gP_M * (1-depMB[t])
          ps[5,t,3,i] <- phiM[i,t] * gP_B * depMB[t] 
          ps[5,t,4,i] <- 0
          ps[5,t,5,i] <- phiM[i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6,i] <- phiM[i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7,i] <- (1 - phiM[i,t])
          
          ps[6,t,1,i] <- 0
          ps[6,t,2,i] <- 0 
          ps[6,t,3,i] <- phiB[i,t] * gP_B 
          ps[6,t,4,i] <- 0
          ps[6,t,5,i] <- 0
          ps[6,t,6,i] <- phiB[i,t] * (1 - gP_B) 
          ps[6,t,7,i] <- (1 - phiB[i,t])
          
          ps[7,t,1,i] <- 0
          ps[7,t,2,i] <- 0
          ps[7,t,3,i] <- 0
          ps[7,t,4,i] <- 0
          ps[7,t,5,i] <- 0
          ps[7,t,6,i] <- 0
          ps[7,t,7,i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1,i] <- pstarW[t]
          po[1,t,2,i] <- 0
          po[1,t,3,i] <- 0
          po[1,t,4,i] <- (1 - pstarW[t])
          
          po[2,t,1,i] <- 0
          po[2,t,2,i] <- pstarM[t]
          po[2,t,3,i] <- 0
          po[2,t,4,i] <- (1 - pstarM[t])
          
          po[3,t,1,i] <- 0
          po[3,t,2,i] <- 0 
          po[3,t,3,i] <- pstarB[t]
          po[3,t,4,i] <- (1 - pstarB[t]) 
          
          po[4,t,1,i] <- 0
          po[4,t,2,i] <- 0
          po[4,t,3,i] <- 0
          po[4,t,4,i] <- 1
          
          po[5,t,1,i] <- 0
          po[5,t,2,i] <- 0
          po[5,t,3,i] <- 0
          po[5,t,4,i] <- 1
            
          po[6,t,1,i] <- 0
          po[6,t,2,i] <- 0
          po[6,t,3,i] <- 0
          po[6,t,4,i] <- 1
          
          po[7,t,1,i] <- 0
          po[7,t,2,i] <- 0
          po[7,t,3,i] <- 0
          po[7,t,4,i] <- 1
        } #i
      } #t



  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()


####-- Winter Covariate Model --####

sink("jags/ms-cjs_W.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  
        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions

            # Winter also varies as individual lm by covariate.
            # No estimates vary by occasion (loop included for future 
            # flexibility)
            
            phiW[i, t] <- mean.phiW[i]
            phiM[i, t] <- mean.phiM
            phiB[i, t] <- mean.phiB
          } #t
        } #i

    # Model for winter covariate

      for(i in 1:n.ind){
        logit(mean.phiW[i]) <- muW[i] + betaW[i]*x[i]
      } #i
    
    muW ~ dnorm(0, 1) # prior for grand mean of logit survival
    betaW ~ dnorm(0, 1) # prior for slope parameter    


    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
        
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #    depMB[i] <- 1
    # }
    
    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    
    # TODO: currently rpoducing impossible beta distributions, fix
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)
    

      # mean.phiW[i] ~ dbeta(5, 0.5)T(0.0001,0.9999)
      
      mean.phiM ~ dbeta(3, 0.5)T(0.0001,0.9999)
      
      mean.phiB ~ dbeta(5, 0.5)T(0.0001,0.9999)


    
    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)
    
    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)

    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)

  
  
  #-- Transition Matrices --#
  
      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1, i] <- phiW[ i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2, i] <- phiW[ i,t] * gPP_M * depWM[t]
          ps[1,t,3, i] <- 0
          ps[1,t,4, i] <- phiW[ i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5, i] <- phiW[ i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6, i] <- 0
          ps[1,t,7, i] <- (1 - phiW[ i,t])
          
          ps[2,t,1, i] <- 0
          ps[2,t,2, i] <- phiM[ i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3, i] <- phiM[ i,t] * gPP_B * depMB[t] 
          ps[2,t,4, i] <- 0
          ps[2,t,5, i] <- phiM[ i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6, i] <- phiM[ i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7, i] <- (1 - phiM[ i,t])
          
          ps[3,t,1, i] <- 0
          ps[3,t,2, i] <- 0 
          ps[3,t,3, i] <- phiB[ i,t] * gPP_B
          ps[3,t,4, i] <- 0 
          ps[3,t,5, i] <- 0
          ps[3,t,6, i] <- phiB[ i,t] * (1 - gPP_B)
          ps[3,t,7, i] <- (1 - phiB[ i,t])
          
          ps[4,t,1, i] <- phiW[ i,t] * gP_W * (1-depWM[t])
          ps[4,t,2, i] <- phiW[ i,t] * gP_M  * depWM[t]
          ps[4,t,3, i] <- 0
          ps[4,t,4, i] <- phiW[ i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5, i] <- phiW[ i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6, i] <- 0
          ps[4,t,7, i] <- (1 - phiW[ i,t]) 
          
          ps[5,t,1, i] <- 0
          ps[5,t,2, i] <- phiM[ i,t] * gP_M * (1-depMB[t])
          ps[5,t,3, i] <- phiM[ i,t] * gP_B * depMB[t] 
          ps[5,t,4, i] <- 0
          ps[5,t,5, i] <- phiM[ i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6, i] <- phiM[ i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7, i] <- (1 - phiM[ i,t])
          
          ps[6,t,1, i] <- 0
          ps[6,t,2, i] <- 0 
          ps[6,t,3, i] <- phiB[ i,t] * gP_B 
          ps[6,t,4, i] <- 0
          ps[6,t,5, i] <- 0
          ps[6,t,6, i] <- phiB[ i,t] * (1 - gP_B) 
          ps[6,t,7, i] <- (1 - phiB[ i,t])
          
          ps[7,t,1, i] <- 0
          ps[7,t,2, i] <- 0
          ps[7,t,3, i] <- 0
          ps[7,t,4, i] <- 0
          ps[7,t,5, i] <- 0
          ps[7,t,6, i] <- 0
          ps[7,t,7, i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1, i] <- pstarW[t]
          po[1,t,2, i] <- 0
          po[1,t,3, i] <- 0
          po[1,t,4, i] <- (1 - pstarW[t])
          
          po[2,t,1, i] <- 0
          po[2,t,2, i] <- pstarM[t]
          po[2,t,3, i] <- 0
          po[2,t,4, i] <- (1 - pstarM[t])
          
          po[3,t,1, i] <- 0
          po[3,t,2, i] <- 0 
          po[3,t,3, i] <- pstarB[t]
          po[3,t,4, i] <- (1 - pstarB[t]) 
          
          po[4,t,1, i] <- 0
          po[4,t,2, i] <- 0
          po[4,t,3, i] <- 0
          po[4,t,4, i] <- 1
          
          po[5,t,1, i] <- 0
          po[5,t,2, i] <- 0
          po[5,t,3, i] <- 0
          po[5,t,4, i] <- 1
            
          po[6,t,1, i] <- 0
          po[6,t,2, i] <- 0
          po[6,t,3, i] <- 0
          po[6,t,4, i] <- 1
          
          po[7,t,1, i] <- 0
          po[7,t,2, i] <- 0
          po[7,t,3, i] <- 0
          po[7,t,4, i] <- 1
        } #i
      } #t

  

  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()


####-- Migration Covariate Model --####

sink("jags/ms-cjs_M.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  
        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions

            # Migration also varies as individual lm by covariate.
            # No estimates vary by occasion (loop included for future 
            # flexibility)
            
            phiW[  i, t] <- mean.phiW
            phiM[  i, t] <- mean.phiM[  i]
            phiB[  i, t] <- mean.phiB
          } #t
        } #i


    # Model for migration covariate


      for(i in 1:n.ind){
        logit(mean.phiM[  i]) <- muM[i] + betaM[i]*x[i]
      } #i
    
    muM ~ dnorm(0, 1) # prior for grand mean of logit survival
    betaM ~ dnorm(0, 1) # prior for slope parameter    


    
    
    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #   depMB[i] <- 1
    # }
    
    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    
    # TODO: currently rpoducing impossible beta distributions, fix
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)
    

      
      mean.phiW ~ dbeta(5, 0.5)T(0.0001,0.9999)
      
      # mean.phiM ~ dbeta(3, 0.5)T(0.0001,0.9999)
      
      mean.phiB ~ dbeta(5, 0.5)T(0.0001,0.9999)


    
    
    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)
    
    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)

    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)

  
  
  #-- Transition Matrices --#
  

      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1, i] <- phiW[ i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2, i] <- phiW[ i,t] * gPP_M * depWM[t]
          ps[1,t,3, i] <- 0
          ps[1,t,4, i] <- phiW[ i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5, i] <- phiW[ i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6, i] <- 0
          ps[1,t,7, i] <- (1 - phiW[ i,t])
          
          ps[2,t,1, i] <- 0
          ps[2,t,2, i] <- phiM[ i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3, i] <- phiM[ i,t] * gPP_B * depMB[t] 
          ps[2,t,4, i] <- 0
          ps[2,t,5, i] <- phiM[ i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6, i] <- phiM[ i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7, i] <- (1 - phiM[ i,t])
          
          ps[3,t,1, i] <- 0
          ps[3,t,2, i] <- 0 
          ps[3,t,3, i] <- phiB[ i,t] * gPP_B
          ps[3,t,4, i] <- 0 
          ps[3,t,5, i] <- 0
          ps[3,t,6, i] <- phiB[ i,t] * (1 - gPP_B)
          ps[3,t,7, i] <- (1 - phiB[ i,t])
          
          ps[4,t,1, i] <- phiW[ i,t] * gP_W * (1-depWM[t])
          ps[4,t,2, i] <- phiW[ i,t] * gP_M  * depWM[t]
          ps[4,t,3, i] <- 0
          ps[4,t,4, i] <- phiW[ i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5, i] <- phiW[ i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6, i] <- 0
          ps[4,t,7, i] <- (1 - phiW[ i,t]) 
          
          ps[5,t,1, i] <- 0
          ps[5,t,2, i] <- phiM[ i,t] * gP_M * (1-depMB[t])
          ps[5,t,3, i] <- phiM[ i,t] * gP_B * depMB[t] 
          ps[5,t,4, i] <- 0
          ps[5,t,5, i] <- phiM[ i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6, i] <- phiM[ i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7, i] <- (1 - phiM[ i,t])
          
          ps[6,t,1, i] <- 0
          ps[6,t,2, i] <- 0 
          ps[6,t,3, i] <- phiB[ i,t] * gP_B 
          ps[6,t,4, i] <- 0
          ps[6,t,5, i] <- 0
          ps[6,t,6, i] <- phiB[ i,t] * (1 - gP_B) 
          ps[6,t,7, i] <- (1 - phiB[ i,t])
          
          ps[7,t,1, i] <- 0
          ps[7,t,2, i] <- 0
          ps[7,t,3, i] <- 0
          ps[7,t,4, i] <- 0
          ps[7,t,5, i] <- 0
          ps[7,t,6, i] <- 0
          ps[7,t,7, i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1, i] <- pstarW[t]
          po[1,t,2, i] <- 0
          po[1,t,3, i] <- 0
          po[1,t,4, i] <- (1 - pstarW[t])
          
          po[2,t,1, i] <- 0
          po[2,t,2, i] <- pstarM[t]
          po[2,t,3, i] <- 0
          po[2,t,4, i] <- (1 - pstarM[t])
          
          po[3,t,1, i] <- 0
          po[3,t,2, i] <- 0 
          po[3,t,3, i] <- pstarB[t]
          po[3,t,4, i] <- (1 - pstarB[t]) 
          
          po[4,t,1, i] <- 0
          po[4,t,2, i] <- 0
          po[4,t,3, i] <- 0
          po[4,t,4, i] <- 1
          
          po[5,t,1, i] <- 0
          po[5,t,2, i] <- 0
          po[5,t,3, i] <- 0
          po[5,t,4, i] <- 1
            
          po[6,t,1, i] <- 0
          po[6,t,2, i] <- 0
          po[6,t,3, i] <- 0
          po[6,t,4, i] <- 1
          
          po[7,t,1, i] <- 0
          po[7,t,2, i] <- 0
          po[7,t,3, i] <- 0
          po[7,t,4, i] <- 1
        } #i
      } #t

  

  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()


####-- Breeding Covariate Model --####

sink("jags/ms-cjs_B.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  
        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions

            # Breeding also varies as individual lm by covariate.
            # No estimates vary by occasion (loop included for future 
            # flexibility)
            
            phiW[  i, t] <- mean.phiW
            phiM[  i, t] <- mean.phiM
            phiB[  i, t] <- mean.phiB[  i]
          } #t
        } #i

    
    # Model for breeding covariate
      for(i in 1:n.ind){
          logit(mean.phiB[ i]) <- muB[i] + betaB[i]*x[i]
        } #i

      muB ~ dnorm(0, 1) # prior for grand mean of logit survival
      betaB ~ dnorm(0, 1) # prior for slope parameter


    
    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #   depMB[i] <- 1
    # }

    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    # TODO: currently producing impossible beta distributions, from the moment
    # matching exercise - fix!
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)
    
      mean.phiW ~ dbeta(5, 0.5)T(0.0001,0.9999)
      
      mean.phiM ~ dbeta(3, 0.5)T(0.0001,0.9999)
      
      # mean.phiB ~ dbeta(5, 0.5)T(0.0001,0.9999)


    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)
        
    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)
    
    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)
  
  
  #-- Transition Matrices --#
  

      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1, i] <- phiW[ i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2, i] <- phiW[ i,t] * gPP_M * depWM[t]
          ps[1,t,3, i] <- 0
          ps[1,t,4, i] <- phiW[ i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5, i] <- phiW[ i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6, i] <- 0
          ps[1,t,7, i] <- (1 - phiW[ i,t])
          
          ps[2,t,1, i] <- 0
          ps[2,t,2, i] <- phiM[ i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3, i] <- phiM[ i,t] * gPP_B * depMB[t] 
          ps[2,t,4, i] <- 0
          ps[2,t,5, i] <- phiM[ i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6, i] <- phiM[ i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7, i] <- (1 - phiM[ i,t])
          
          ps[3,t,1, i] <- 0
          ps[3,t,2, i] <- 0 
          ps[3,t,3, i] <- phiB[ i,t] * gPP_B
          ps[3,t,4, i] <- 0 
          ps[3,t,5, i] <- 0
          ps[3,t,6, i] <- phiB[ i,t] * (1 - gPP_B)
          ps[3,t,7, i] <- (1 - phiB[ i,t])
          
          ps[4,t,1, i] <- phiW[ i,t] * gP_W * (1-depWM[t])
          ps[4,t,2, i] <- phiW[ i,t] * gP_M  * depWM[t]
          ps[4,t,3, i] <- 0
          ps[4,t,4, i] <- phiW[ i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5, i] <- phiW[ i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6, i] <- 0
          ps[4,t,7, i] <- (1 - phiW[ i,t]) 
          
          ps[5,t,1, i] <- 0
          ps[5,t,2, i] <- phiM[ i,t] * gP_M * (1-depMB[t])
          ps[5,t,3, i] <- phiM[ i,t] * gP_B * depMB[t] 
          ps[5,t,4, i] <- 0
          ps[5,t,5, i] <- phiM[ i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6, i] <- phiM[ i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7, i] <- (1 - phiM[ i,t])
          
          ps[6,t,1, i] <- 0
          ps[6,t,2, i] <- 0 
          ps[6,t,3, i] <- phiB[ i,t] * gP_B 
          ps[6,t,4, i] <- 0
          ps[6,t,5, i] <- 0
          ps[6,t,6, i] <- phiB[ i,t] * (1 - gP_B) 
          ps[6,t,7, i] <- (1 - phiB[ i,t])
          
          ps[7,t,1, i] <- 0
          ps[7,t,2, i] <- 0
          ps[7,t,3, i] <- 0
          ps[7,t,4, i] <- 0
          ps[7,t,5, i] <- 0
          ps[7,t,6, i] <- 0
          ps[7,t,7, i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1, i] <- pstarW[t]
          po[1,t,2, i] <- 0
          po[1,t,3, i] <- 0
          po[1,t,4, i] <- (1 - pstarW[t])
          
          po[2,t,1, i] <- 0
          po[2,t,2, i] <- pstarM[t]
          po[2,t,3, i] <- 0
          po[2,t,4, i] <- (1 - pstarM[t])
          
          po[3,t,1, i] <- 0
          po[3,t,2, i] <- 0 
          po[3,t,3, i] <- pstarB[t]
          po[3,t,4, i] <- (1 - pstarB[t]) 
          
          po[4,t,1, i] <- 0
          po[4,t,2, i] <- 0
          po[4,t,3, i] <- 0
          po[4,t,4, i] <- 1
          
          po[5,t,1, i] <- 0
          po[5,t,2, i] <- 0
          po[5,t,3, i] <- 0
          po[5,t,4, i] <- 1
            
          po[6,t,1, i] <- 0
          po[6,t,2, i] <- 0
          po[6,t,3, i] <- 0
          po[6,t,4, i] <- 1
          
          po[7,t,1, i] <- 0
          po[7,t,2, i] <- 0
          po[7,t,3, i] <- 0
          po[7,t,4, i] <- 1
        } #i
      } #t

  

  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()

####-- Winter & Migration Covariate Model --####

sink("jags/ms-cjs_WM.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  

        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions

            # Winter & Migration also vary as individual lm by covariate.
            # No estimates vary by occasion (loop included for future 
            # flexibility)
            
            phiW[  i, t] <- mean.phiW[  i]
            phiM[  i, t] <- mean.phiM[  i]
            phiB[  i, t] <- mean.phiB
          } #t
        } #i


    # Model for winter covariate


      for(i in 1:n.ind){
        logit(mean.phiW[  i]) <- muW[i] + betaW[i]*x[i]
      } #i
    
    muW ~ dnorm(0, 1) # prior for grand mean of logit survival
    betaW ~ dnorm(0, 1) # prior for slope parameter    

    
    # Model for migration covariate


      for(i in 1:n.ind){
        logit(mean.phiM[i]) <- muM[i] + betaM[i]*x[i]
      } #i
    
    muM ~ dnorm(0, 1) # prior for grand mean of logit survival
    betaM ~ dnorm(0, 1) # prior for slope parameter    


    
    
    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
        
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #   depMB[i] <- 1
    # }

    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    
    # TODO: currently rpoducing impossible beta distributions, fix
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)

      
      # mean.phiW ~ dbeta(5, 0.5)T(0.0001,0.9999)
      
      # mean.phiM ~ dbeta(3, 0.5)T(0.0001,0.9999)
      
      mean.phiB ~ dbeta(5, 0.5)T(0.0001,0.9999)

    
    
    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)
          
    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)
    
    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)

  
  #-- Transition Matrices --#

      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1, i] <- phiW[ i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2, i] <- phiW[ i,t] * gPP_M * depWM[t]
          ps[1,t,3, i] <- 0
          ps[1,t,4, i] <- phiW[ i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5, i] <- phiW[ i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6, i] <- 0
          ps[1,t,7, i] <- (1 - phiW[ i,t])
          
          ps[2,t,1, i] <- 0
          ps[2,t,2, i] <- phiM[ i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3, i] <- phiM[ i,t] * gPP_B * depMB[t] 
          ps[2,t,4, i] <- 0
          ps[2,t,5, i] <- phiM[ i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6, i] <- phiM[ i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7, i] <- (1 - phiM[ i,t])
          
          ps[3,t,1, i] <- 0
          ps[3,t,2, i] <- 0 
          ps[3,t,3, i] <- phiB[ i,t] * gPP_B
          ps[3,t,4, i] <- 0 
          ps[3,t,5, i] <- 0
          ps[3,t,6, i] <- phiB[ i,t] * (1 - gPP_B)
          ps[3,t,7, i] <- (1 - phiB[ i,t])
          
          ps[4,t,1, i] <- phiW[ i,t] * gP_W * (1-depWM[t])
          ps[4,t,2, i] <- phiW[ i,t] * gP_M  * depWM[t]
          ps[4,t,3, i] <- 0
          ps[4,t,4, i] <- phiW[ i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5, i] <- phiW[ i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6, i] <- 0
          ps[4,t,7, i] <- (1 - phiW[ i,t]) 
          
          ps[5,t,1, i] <- 0
          ps[5,t,2, i] <- phiM[ i,t] * gP_M * (1-depMB[t])
          ps[5,t,3, i] <- phiM[ i,t] * gP_B * depMB[t] 
          ps[5,t,4, i] <- 0
          ps[5,t,5, i] <- phiM[ i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6, i] <- phiM[ i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7, i] <- (1 - phiM[ i,t])
          
          ps[6,t,1, i] <- 0
          ps[6,t,2, i] <- 0 
          ps[6,t,3, i] <- phiB[ i,t] * gP_B 
          ps[6,t,4, i] <- 0
          ps[6,t,5, i] <- 0
          ps[6,t,6, i] <- phiB[ i,t] * (1 - gP_B) 
          ps[6,t,7, i] <- (1 - phiB[ i,t])
          
          ps[7,t,1, i] <- 0
          ps[7,t,2, i] <- 0
          ps[7,t,3, i] <- 0
          ps[7,t,4, i] <- 0
          ps[7,t,5, i] <- 0
          ps[7,t,6, i] <- 0
          ps[7,t,7, i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1, i] <- pstarW[t]
          po[1,t,2, i] <- 0
          po[1,t,3, i] <- 0
          po[1,t,4, i] <- (1 - pstarW[t])
          
          po[2,t,1, i] <- 0
          po[2,t,2, i] <- pstarM[t]
          po[2,t,3, i] <- 0
          po[2,t,4, i] <- (1 - pstarM[t])
          
          po[3,t,1, i] <- 0
          po[3,t,2, i] <- 0 
          po[3,t,3, i] <- pstarB[t]
          po[3,t,4, i] <- (1 - pstarB[t]) 
          
          po[4,t,1, i] <- 0
          po[4,t,2, i] <- 0
          po[4,t,3, i] <- 0
          po[4,t,4, i] <- 1
          
          po[5,t,1, i] <- 0
          po[5,t,2, i] <- 0
          po[5,t,3, i] <- 0
          po[5,t,4, i] <- 1
            
          po[6,t,1, i] <- 0
          po[6,t,2, i] <- 0
          po[6,t,3, i] <- 0
          po[6,t,4, i] <- 1
          
          po[7,t,1, i] <- 0
          po[7,t,2, i] <- 0
          po[7,t,3, i] <- 0
          po[7,t,4, i] <- 1
        } #i
      } #t

  

  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()

####-- Winter & Breeding Covariate Model --####

sink("jags/ms-cjs_WB.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  

        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions

            # Winter & Breeding also vary as individual lm by covariate.
            # No estimates vary by occasion (loop included for future 
            # flexibility)
            
            phiW[  i, t] <- mean.phiW[  i]
            phiM[  i, t] <- mean.phiM
            phiB[  i, t] <- mean.phiB[  i]
          } #t
        } #i

    
    # Model for winter covariate


      for(i in 1:n.ind){
        logit(mean.phiW[i]) <- muW[i] + betaW[i]*x[i]
      } #i
    
    muW ~ dnorm(0, 1) # prior for grand mean of logit survival
    betaW ~ dnorm(0, 1) # prior for slope parameter    


    
    # Model for breeding covariate


        for(i in 1:n.ind){
          logit(mean.phiB[i]) <- muB[i] + betaB[i]*x[i]
        } #i

      muB ~ dnorm(0, 1) # prior for grand mean of logit survival
      betaB ~ dnorm(0, 1) # prior for slope parameter


    
    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
        
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #   depMB[i] <- 1
    # }

    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    
    # TODO: currently rpoducing impossible beta distributions, fix
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)
    

      
      # mean.phiW ~ dbeta(5, 0.5)T(0.0001,0.9999)
      
      mean.phiM ~ dbeta(3, 0.5)T(0.0001,0.9999)
      
      # mean.phiB ~ dbeta(5, 0.5)T(0.0001,0.9999)


    
    
    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)
            
    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)
    
    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)

  
  #-- Transition Matrices --#

      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1, i] <- phiW[ i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2, i] <- phiW[ i,t] * gPP_M * depWM[t]
          ps[1,t,3, i] <- 0
          ps[1,t,4, i] <- phiW[ i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5, i] <- phiW[ i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6, i] <- 0
          ps[1,t,7, i] <- (1 - phiW[ i,t])
          
          ps[2,t,1, i] <- 0
          ps[2,t,2, i] <- phiM[ i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3, i] <- phiM[ i,t] * gPP_B * depMB[t] 
          ps[2,t,4, i] <- 0
          ps[2,t,5, i] <- phiM[ i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6, i] <- phiM[ i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7, i] <- (1 - phiM[ i,t])
          
          ps[3,t,1, i] <- 0
          ps[3,t,2, i] <- 0 
          ps[3,t,3, i] <- phiB[ i,t] * gPP_B
          ps[3,t,4, i] <- 0 
          ps[3,t,5, i] <- 0
          ps[3,t,6, i] <- phiB[ i,t] * (1 - gPP_B)
          ps[3,t,7, i] <- (1 - phiB[ i,t])
          
          ps[4,t,1, i] <- phiW[ i,t] * gP_W * (1-depWM[t])
          ps[4,t,2, i] <- phiW[ i,t] * gP_M  * depWM[t]
          ps[4,t,3, i] <- 0
          ps[4,t,4, i] <- phiW[ i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5, i] <- phiW[ i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6, i] <- 0
          ps[4,t,7, i] <- (1 - phiW[ i,t]) 
          
          ps[5,t,1, i] <- 0
          ps[5,t,2, i] <- phiM[ i,t] * gP_M * (1-depMB[t])
          ps[5,t,3, i] <- phiM[ i,t] * gP_B * depMB[t] 
          ps[5,t,4, i] <- 0
          ps[5,t,5, i] <- phiM[ i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6, i] <- phiM[ i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7, i] <- (1 - phiM[ i,t])
          
          ps[6,t,1, i] <- 0
          ps[6,t,2, i] <- 0 
          ps[6,t,3, i] <- phiB[ i,t] * gP_B 
          ps[6,t,4, i] <- 0
          ps[6,t,5, i] <- 0
          ps[6,t,6, i] <- phiB[ i,t] * (1 - gP_B) 
          ps[6,t,7, i] <- (1 - phiB[ i,t])
          
          ps[7,t,1, i] <- 0
          ps[7,t,2, i] <- 0
          ps[7,t,3, i] <- 0
          ps[7,t,4, i] <- 0
          ps[7,t,5, i] <- 0
          ps[7,t,6, i] <- 0
          ps[7,t,7, i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1, i] <- pstarW[t]
          po[1,t,2, i] <- 0
          po[1,t,3, i] <- 0
          po[1,t,4, i] <- (1 - pstarW[t])
          
          po[2,t,1, i] <- 0
          po[2,t,2, i] <- pstarM[t]
          po[2,t,3, i] <- 0
          po[2,t,4, i] <- (1 - pstarM[t])
          
          po[3,t,1, i] <- 0
          po[3,t,2, i] <- 0 
          po[3,t,3, i] <- pstarB[t]
          po[3,t,4, i] <- (1 - pstarB[t]) 
          
          po[4,t,1, i] <- 0
          po[4,t,2, i] <- 0
          po[4,t,3, i] <- 0
          po[4,t,4, i] <- 1
          
          po[5,t,1, i] <- 0
          po[5,t,2, i] <- 0
          po[5,t,3, i] <- 0
          po[5,t,4, i] <- 1
            
          po[6,t,1, i] <- 0
          po[6,t,2, i] <- 0
          po[6,t,3, i] <- 0
          po[6,t,4, i] <- 1
          
          po[7,t,1, i] <- 0
          po[7,t,2, i] <- 0
          po[7,t,3, i] <- 0
          po[7,t,4, i] <- 1
        } #i
      } #t

  

  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()


####-- Migration & Breeding Covariate Model --####

sink("jags/ms-cjs_MB.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  

        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions

            # Migration & Breeding also vary as individual lm by covariate.
            # No estimates vary by occasion (loop included for future 
            # flexibility)
            
            phiW[  i, t] <- mean.phiW
            phiM[  i, t] <- mean.phiM[  i]
            phiB[  i, t] <- mean.phiB[  i]
          } #t
        } #i

    
    # Model for migration covariate


      for(i in 1:n.ind){
        logit(mean.phiM[i]) <- muM[i] + betaM[i]*x[i]
      } #i
    
    muM ~ dnorm(0, 1) # prior for grand mean of logit survival
    betaM ~ dnorm(0, 1) # prior for slope parameter    

    
    # Model for breeding covariate


        for(i in 1:n.ind){
          logit(mean.phiB[i]) <- muB[i] + betaB[i]*x[i]
        } #i

      muB ~ dnorm(0, 1) # prior for grand mean of logit survival
      betaB ~ dnorm(0, 1) # prior for slope parameter

    
    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #   depMB[i] <- 1
    # }

    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    
    # TODO: currently rpoducing impossible beta distributions, fix
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)
    

      
      mean.phiW ~ dbeta(5, 0.5)T(0.0001,0.9999)
      
      # mean.phiM ~ dbeta(3, 0.5)T(0.0001,0.9999)
      
      # mean.phiB ~ dbeta(5, 0.5)T(0.0001,0.9999)


    
    
    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)
        
    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)
    
    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)
  
  
  #-- Transition Matrices --#
  

      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1, i] <- phiW[ i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2, i] <- phiW[ i,t] * gPP_M * depWM[t]
          ps[1,t,3, i] <- 0
          ps[1,t,4, i] <- phiW[ i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5, i] <- phiW[ i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6, i] <- 0
          ps[1,t,7, i] <- (1 - phiW[ i,t])
          
          ps[2,t,1, i] <- 0
          ps[2,t,2, i] <- phiM[ i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3, i] <- phiM[ i,t] * gPP_B * depMB[t] 
          ps[2,t,4, i] <- 0
          ps[2,t,5, i] <- phiM[ i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6, i] <- phiM[ i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7, i] <- (1 - phiM[ i,t])
          
          ps[3,t,1, i] <- 0
          ps[3,t,2, i] <- 0 
          ps[3,t,3, i] <- phiB[ i,t] * gPP_B
          ps[3,t,4, i] <- 0 
          ps[3,t,5, i] <- 0
          ps[3,t,6, i] <- phiB[ i,t] * (1 - gPP_B)
          ps[3,t,7, i] <- (1 - phiB[ i,t])
          
          ps[4,t,1, i] <- phiW[ i,t] * gP_W * (1-depWM[t])
          ps[4,t,2, i] <- phiW[ i,t] * gP_M  * depWM[t]
          ps[4,t,3, i] <- 0
          ps[4,t,4, i] <- phiW[ i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5, i] <- phiW[ i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6, i] <- 0
          ps[4,t,7, i] <- (1 - phiW[ i,t]) 
          
          ps[5,t,1, i] <- 0
          ps[5,t,2, i] <- phiM[ i,t] * gP_M * (1-depMB[t])
          ps[5,t,3, i] <- phiM[ i,t] * gP_B * depMB[t] 
          ps[5,t,4, i] <- 0
          ps[5,t,5, i] <- phiM[ i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6, i] <- phiM[ i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7, i] <- (1 - phiM[ i,t])
          
          ps[6,t,1, i] <- 0
          ps[6,t,2, i] <- 0 
          ps[6,t,3, i] <- phiB[ i,t] * gP_B 
          ps[6,t,4, i] <- 0
          ps[6,t,5, i] <- 0
          ps[6,t,6, i] <- phiB[ i,t] * (1 - gP_B) 
          ps[6,t,7, i] <- (1 - phiB[ i,t])
          
          ps[7,t,1, i] <- 0
          ps[7,t,2, i] <- 0
          ps[7,t,3, i] <- 0
          ps[7,t,4, i] <- 0
          ps[7,t,5, i] <- 0
          ps[7,t,6, i] <- 0
          ps[7,t,7, i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1, i] <- pstarW[t]
          po[1,t,2, i] <- 0
          po[1,t,3, i] <- 0
          po[1,t,4, i] <- (1 - pstarW[t])
          
          po[2,t,1, i] <- 0
          po[2,t,2, i] <- pstarM[t]
          po[2,t,3, i] <- 0
          po[2,t,4, i] <- (1 - pstarM[t])
          
          po[3,t,1, i] <- 0
          po[3,t,2, i] <- 0 
          po[3,t,3, i] <- pstarB[t]
          po[3,t,4, i] <- (1 - pstarB[t]) 
          
          po[4,t,1, i] <- 0
          po[4,t,2, i] <- 0
          po[4,t,3, i] <- 0
          po[4,t,4, i] <- 1
          
          po[5,t,1, i] <- 0
          po[5,t,2, i] <- 0
          po[5,t,3, i] <- 0
          po[5,t,4, i] <- 1
            
          po[6,t,1, i] <- 0
          po[6,t,2, i] <- 0
          po[6,t,3, i] <- 0
          po[6,t,4, i] <- 1
          
          po[7,t,1, i] <- 0
          po[7,t,2, i] <- 0
          po[7,t,3, i] <- 0
          po[7,t,4, i] <- 1
        } #i
      } #t

  

  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()


####-- Winter, Migration, & Breeding Covariate Model --####

sink("jags/ms-cjs_WMB.jags")
cat("
model{
  
  #-- Priors --#
  
    #- Probability of Survival (phi)
  

        for(i in 1:n.ind){
          for(t in 1:(n.prim-1)){ # loop through occasions

            # Winter, Migration, & Breeding also vary as individual lm by covariate.
            # No estimates vary by occasion (loop included for future 
            # flexibility)
            
            phiW[  i, t] <- mean.phiW[  i]
            phiM[  i, t] <- mean.phiM[  i]
            phiB[  i, t] <- mean.phiB[  i]
          } #t
        } #i

    
    # Model for winter covariate


      for(i in 1:n.ind){
        logit(mean.phiW[i]) <- muW + betaW*x[i]
      } #i
    
    muW~ dnorm(0, 1) # prior for grand mean of logit survival
    betaW ~ dnorm(0, 1) # prior for slope parameter    


    
    # Model for migration covariate


      for(i in 1:n.ind){
        logit(mean.phiM[i]) <- muM + betaM*x[i]
      } #i
    
    muM ~ dnorm(0, 1) # prior for grand mean of logit survival
    betaM ~ dnorm(0, 1) # prior for slope parameter    

    
    # Model for breeding covariate


        for(i in 1:n.ind){
          logit(mean.phiB[i]) <- muB + betaB*x[i]
        } #i

      muB ~ dnorm(0, 1) # prior for grand mean of logit survival
      betaB ~ dnorm(0, 1) # prior for slope parameter

    
    #- Probability of Detection (p) for secondaries
    
    for (t in 1:n.prim){ # loop thorugh primary occasions
      for (j in 1:max(n.sec[1:n.prim])){ # loop through secondar occasions
      
        # constant probability of detection
        
        pW[t,j] <- mean.pW
        pM[t,j] <- mean.pM
        pB[t,j] <- mean.pB
      } #j
    } #t
    
    
    #- Bernoulli trial within secondaries
    
    for (t in 1:n.prim){
      for (j in 1:n.sec[t]){
        yesW[t,j] ~ dbin(pW[t,j], totalW[t,j])
        yesM[t,j] ~ dbin(pM[t,j], totalM[t,j])
        yesB[t,j] ~ dbin(pB[t,j], totalB[t,j])
      }
    }
    
    
    #- Primary occasion pooled detection probability
    
    for (t in 1:n.prim){
      pstarW[t] <- 1 - prod(1 - pW[t,1:n.sec[t]])
      pstarM[t] <- 1 - prod(1 - pM[t,1:n.sec[t]])
      pstarB[t] <- 1 - prod(1 - pB[t,1:n.sec[t]])
    }
  
    
    #- Seasonal Transition Probabilities
  
    # Departure probabilities (based on method in Pledger et al. 2009)
    depWM[1] <- 0 #set initial departure prob to 0
    
    # No transition to breeding until first obs of any ind in migration
    for(i in 1:5){
      depMB[i] <- 0 
    }
    
    # # Enforce migration completed by weeks 16-17 (~June 12)
    # for(i in (n.prim-1):n.prim){
    #   depMB[i] <- 1
    # }

    # Set departure probs to be estimated as a Weibull CDF
    
    # wint -> mig
    for(t in 2:n.prim) {
      depWM[t] <- 1-exp(-(t/gWM)^kWM + ((t-1)/gWM)^kWM)
    }
    
    # mig -> breed
    for(t in 6:n.prim) {
      depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
    }
  
  
  #-- Hyperparameters --#
  
    #- Phi
    
    # (derived from moment matching based on estimates in Rockwell et al. 2017)
    
    # TODO: currently rpoducing impossible beta distributions, fix
    # mean.phiM ~ dbeta(1.736539, 0.05131568)T(0.0001,0.9999)
    # mean.phiW ~ dbeta(1.598064, 0.008418654)T(0.0001,0.9999)
    # mean.phiB ~ dbeta(3.175514, 0.02714958)T(0.0001,0.9999)

    #- P
    
    mean.pW ~ dbeta(1, 1)
    mean.pM ~ dbeta(1, 1)
    mean.pB ~ dbeta(1, 1)
    
    
    #- Markovian availability
    
    # Probability of availability given prevously unobserved
    gP_W ~ dbeta(1, 1)
    gP_M ~ dbeta(1, 1)
    gP_B ~ dbeta(1, 1)
    
    # Probability of availability given prevously observed
    gPP_W ~ dbeta(1, 1)
    gPP_M ~ dbeta(1, 1)
    gPP_B ~ dbeta(1, 1)

    # Hyperparams for the Weibull params
    gWM ~ dunif(4, 12)
    gMB ~ dunif(8, 16)
    
    kWM ~ dgamma(10, 1)
    kMB ~ dgamma(10, 1)
      
  
   #-- Transition Matrices --#

      for(t in 1:(n.prim-1)){
        for(i in 1:n.ind){
          #- State Transitons
          
          ps[1,t,1, i] <- phiW[ i,t] * gPP_W * (1-depWM[t])
          ps[1,t,2, i] <- phiW[ i,t] * gPP_M * depWM[t]
          ps[1,t,3, i] <- 0
          ps[1,t,4, i] <- phiW[ i,t] * (1- gPP_W) * (1-depWM[t])
          ps[1,t,5, i] <- phiW[ i,t] * (1- gPP_M) * (depWM[t])
          ps[1,t,6, i] <- 0
          ps[1,t,7, i] <- (1 - phiW[ i,t])
          
          ps[2,t,1, i] <- 0
          ps[2,t,2, i] <- phiM[ i,t] * gPP_M * (1-depMB[t])
          ps[2,t,3, i] <- phiM[ i,t] * gPP_B * depMB[t] 
          ps[2,t,4, i] <- 0
          ps[2,t,5, i] <- phiM[ i,t] * (1 - gPP_M) * (1-depMB[t]) 
          ps[2,t,6, i] <- phiM[ i,t] * (1 - gPP_B) * depMB[t]
          ps[2,t,7, i] <- (1 - phiM[ i,t])
          
          ps[3,t,1, i] <- 0
          ps[3,t,2, i] <- 0 
          ps[3,t,3, i] <- phiB[ i,t] * gPP_B
          ps[3,t,4, i] <- 0 
          ps[3,t,5, i] <- 0
          ps[3,t,6, i] <- phiB[ i,t] * (1 - gPP_B)
          ps[3,t,7, i] <- (1 - phiB[ i,t])
          
          ps[4,t,1, i] <- phiW[ i,t] * gP_W * (1-depWM[t])
          ps[4,t,2, i] <- phiW[ i,t] * gP_M  * depWM[t]
          ps[4,t,3, i] <- 0
          ps[4,t,4, i] <- phiW[ i,t] * (1- gP_W) * (1-depWM[t])
          ps[4,t,5, i] <- phiW[ i,t] * (1- gP_M) * depWM[t]
          ps[4,t,6, i] <- 0
          ps[4,t,7, i] <- (1 - phiW[ i,t]) 
          
          ps[5,t,1, i] <- 0
          ps[5,t,2, i] <- phiM[ i,t] * gP_M * (1-depMB[t])
          ps[5,t,3, i] <- phiM[ i,t] * gP_B * depMB[t] 
          ps[5,t,4, i] <- 0
          ps[5,t,5, i] <- phiM[ i,t] * (1 - gP_M) * (1-depMB[t])
          ps[5,t,6, i] <- phiM[ i,t] * (1 - gP_B) * depMB[t] 
          ps[5,t,7, i] <- (1 - phiM[ i,t])
          
          ps[6,t,1, i] <- 0
          ps[6,t,2, i] <- 0 
          ps[6,t,3, i] <- phiB[ i,t] * gP_B 
          ps[6,t,4, i] <- 0
          ps[6,t,5, i] <- 0
          ps[6,t,6, i] <- phiB[ i,t] * (1 - gP_B) 
          ps[6,t,7, i] <- (1 - phiB[ i,t])
          
          ps[7,t,1, i] <- 0
          ps[7,t,2, i] <- 0
          ps[7,t,3, i] <- 0
          ps[7,t,4, i] <- 0
          ps[7,t,5, i] <- 0
          ps[7,t,6, i] <- 0
          ps[7,t,7, i] <- 1
          
          
          #- Observation Model
          
          po[1,t,1, i] <- pstarW[t]
          po[1,t,2, i] <- 0
          po[1,t,3, i] <- 0
          po[1,t,4, i] <- (1 - pstarW[t])
          
          po[2,t,1, i] <- 0
          po[2,t,2, i] <- pstarM[t]
          po[2,t,3, i] <- 0
          po[2,t,4, i] <- (1 - pstarM[t])
          
          po[3,t,1, i] <- 0
          po[3,t,2, i] <- 0 
          po[3,t,3, i] <- pstarB[t]
          po[3,t,4, i] <- (1 - pstarB[t]) 
          
          po[4,t,1, i] <- 0
          po[4,t,2, i] <- 0
          po[4,t,3, i] <- 0
          po[4,t,4, i] <- 1
          
          po[5,t,1, i] <- 0
          po[5,t,2, i] <- 0
          po[5,t,3, i] <- 0
          po[5,t,4, i] <- 1
            
          po[6,t,1, i] <- 0
          po[6,t,2, i] <- 0
          po[6,t,3, i] <- 0
          po[6,t,4, i] <- 1
          
          po[7,t,1, i] <- 0
          po[7,t,2, i] <- 0
          po[7,t,3, i] <- 0
          po[7,t,4, i] <- 1
        } #i
      } #t

  

  ##-- Likelihood --##
  
  for(i in 1:n.ind){
    
    # Known to be alive and in winter state on first occasion
    # z[i, first[i]] <- ch[i, first[i]]
    
    for(t in (first[i] + 1):n.prim){
      
      # State model
      z[i, t] ~ dcat(ps[z[i, t-1], t-1,  , i])
      
      # Observation model
      ch[i, t] ~ dcat(po[z[i, t], t-1,  , i])
    } #t
  } #i
  
} #model
    ", fill = TRUE)
sink()

##---- Finalize Script ----##
message("Script complete.")
