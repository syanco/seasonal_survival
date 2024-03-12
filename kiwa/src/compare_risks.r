# Script to generate relative risks and odds ratios (with uncertainty) for two
# survival estimates.


#---- INITS ----#

# Libraries
library(tidyverse)
library(glue)
library(emmeans)
library(rjags)

# Functions
# function anticipates a data frame with 2 rows and the following columns:
# 1 name: the name of each survival estimate to be compared
# 2 est: the survival estimate
# 3 lower: the lower uncertainty bound (eg CI) for the survival estimate
# 4 upper: the upper uncertainty bound for the survival estimate
# 5 dur: exposure duration in units of the survival estimate (e.g., if estimated 
#       survival is daily and the period of interst is one week, the dur = 7)

compareSurv <- function(data, RR = T, OR = T, with_dur = T){
  if(nrow(data) != 2){
    message("Data must be dataframe with two rows.")
    return(NULL)
  }
  
  # unpack estimates
  if(with_dur == T){
    pdeath1_est <- 1-(data$est[1])^data$dur[1]
    pdeath1_low <- 1-(data$lower[1])^data$dur[1]
    pdeath1_hi <- 1-(data$upper[1])^data$dur[1]
    pdeath2_est <- 1-(data$est[2])^data$dur[2]
    pdeath2_low <- 1-(data$lower[2])^data$dur[2]
    pdeath2_hi <- 1-(data$upper[2])^data$dur[2]
  } else{
    pdeath1_est <- 1-data$est[1]
    pdeath1_low <- 1-data$lower[1]
    pdeath1_hi <- 1-data$upper[1]
    pdeath2_est <- 1-data$est[2]
    pdeath2_low <- 1-data$lower[2]
    pdeath2_hi <- 1-data$upper[2]
  }
  # estimated comparison
  (RR_est <- pdeath1_est/pdeath2_est)
  
  (odds1_est <- pdeath1_est/(1-pdeath1_est))
  (odds2_est <- pdeath2_est/(1-pdeath2_est))
  
  (OR_est <- odds1_est/odds2_est) 
  
  
  # low v high comparison
  (RR_lh <- pdeath1_low/pdeath2_hi)
  
  (odds1_l <- pdeath1_low/(1-pdeath1_low))
  (odds2_h <- pdeath2_hi/(1-pdeath2_hi))
  
  (OR_lh <- odds1_l/odds2_h)
  
  
  # high v low comparison
  (RR_hl <- pdeath1_hi/pdeath2_low)
  
  (odds1_h <- pdeath1_hi/(1-pdeath1_hi))
  (odds2_l <- pdeath2_low/(1-pdeath2_low))
  
  (OR_hl <- odds1_h/odds2_l)
  
  # collate output
  out <- data.frame(comparison = glue("{data$name[1]}/{data$name[2]}"),
                    estimated_RR = RR_est,
                    LvH_RR = RR_lh,
                    HvL_RR = RR_hl,
                    estimated_OR = OR_est,
                    LvH_OR = OR_lh,
                    HvL_OR = OR_hl)
  
  print(out)
  return(out)
}


getConditionalSurv <- function(model, covar, intercept, slope, foc_season){
  m <- c()
  quant_l <- c()
  quant_h <- c()
  # age_l <- c()
  post <- plogis(intercept+slope*covar)
  m <- mean(post)
  quant_l <- quantile(post, probs = c(0.025))
  quant_h <- quantile(post, probs = c(0.975))
  # age_l <- age_class
  season <- foc_season
  
  cond  <- data.frame(estimate = m,
                      ci_l = quant_l,
                      ci_h = quant_h,
                      # age = age_l,
                      season = season)
  return(cond)
}



#---- LOAD DATA ----#

# get list of model rdata files
mods <- list.files("output", full.names = T, pattern = "noage")
mods_nums <- c(1,2)

# stupid for loop I had to make b/c I did a bad job naming model objects upon creation...
for(i in 1:length(mods)){
  load(mods[i])
  assign(glue("mod{mods_nums[i]}"), mod)
}

#define focal mod
# chosen in other results script
focal_mod <- mod2

# Load EVI Data
load("output/kiwa_dat_JAGS.rdata")
raws <- data.frame(evi = evi)

min_evi <- min(na.omit(evi))
max_evi <- max(na.omit(evi))
mean_evi <- mean(na.omit(evi))
eviq <- quantile(evi, na.rm = T, probs = c(0.05, 0.95))



#---- CALCULATE CONTRASTS ----#

#-- Compare Seasons at mean EVI

#- Get conditional estimates at mean EVI
# Breeding
cond_b <- getConditionalSurv(model = focal_mod, covar = mean_evi, 
                   intercept = focal_mod$sims.list$muB,
                   slope = focal_mod$sims.list$betaB,
                   foc_season = "breeding")

# Winter
cond_w  <- getConditionalSurv(model = focal_mod, covar = mean_evi, 
                              intercept = focal_mod$sims.list$muW,
                              slope = focal_mod$sims.list$betaW,
                              foc_season = "winter")

# Migration
cond_m  <- getConditionalSurv(model = focal_mod, covar = mean_evi, 
                              intercept = focal_mod$sims.list$muM,
                              slope = focal_mod$sims.list$betaM,
                              foc_season = "migration")

# Bind together
conds <- rbind(cond_b, cond_w, cond_m) %>% 
  mutate(dur = case_when(season == "breeding" ~ (130/7),
                         season == "migration" ~ (17/7),
                         season == "winter" ~ (189/7)))

#- Compare across seasons

# Breeding v. Migration
data_bm <- conds %>% 
  filter(season == "breeding" | season == "migration") %>% 
  rename(name = season, 
         est = estimate,
         upper = ci_h,
         lower = ci_l,
         dur = dur) 

bm <- compareSurv(data=data_bm)

# Breeding v. Winter
data_bw <- conds %>% 
  filter(season == "breeding" | season == "winter") %>% 
  rename(name = season, 
         est = estimate,
         upper = ci_h,
         lower = ci_l,
         dur = dur) 

bw <- compareSurv(data=data_bw)

# Winter v. Migration
data_wm <- conds %>% 
  filter(season == "winter" | season == "migration") %>% 
  rename(name = season, 
         est = estimate,
         upper = ci_h,
         lower = ci_l,
         dur = dur) 

wm <- compareSurv(data=data_wm)

# Bind results and write out
comp_seasons <- rbind(bm, bw, wm)
write.csv(comp_seasons, file = "output/comp_season_risk.csv")

#---- COMPARE ACROSS EVI (W/I SEASON) ----#

#-- winter

# get Conditional estimates at low (5% quantile) and high (95% quantile) EVI values

# low evi
cond_wl <- getConditionalSurv(model = focal_mod, covar = eviq[1], 
                   intercept = focal_mod$sims.list$muW,
                   slope = focal_mod$sims.list$betaW,
                   foc_season = "winter_low_evi")
# high evi
cond_wh<- getConditionalSurv(model = focal_mod, covar = eviq[2], 
                              intercept = focal_mod$sims.list$muW,
                              slope = focal_mod$sims.list$betaW,
                              foc_season = "winter_high_evi")

# Bind results
data_wint <- rbind(cond_wl, cond_wh) %>% 
  rename(name = season, 
         est = estimate,
         upper = ci_h,
         lower = ci_l) %>% 
         mutate(dur = (189/7)) 

wint_evi <- compareSurv(data=data_wint)


#-- breeding

# get Conditional estimates at low (5% quantile) and high (95% quantile) EVI values

# low evi
cond_bl <- getConditionalSurv(model = focal_mod, covar = eviq[1], 
                              intercept = focal_mod$sims.list$muB,
                              slope = focal_mod$sims.list$betaB,
                              foc_season = "breeding_low_evi")
# high evi
cond_bh<- getConditionalSurv(model = focal_mod, covar = eviq[2], 
                             intercept = focal_mod$sims.list$muB,
                             slope = focal_mod$sims.list$betaB,
                             foc_season = "breeding_high_evi")

# Bind results
data_breed <- rbind(cond_bl, cond_bh) %>% 
  rename(name = season, 
         est = estimate,
         upper = ci_h,
         lower = ci_l) %>% 
  mutate(dur = (130/7)) 

breed_evi <- compareSurv(data=data_breed)


#-- Migration

# get Conditional estimates at low (5% quantile) and high (95% quantile) EVI values

# low evi
cond_ml <- getConditionalSurv(model = focal_mod, covar = eviq[1], 
                              intercept = focal_mod$sims.list$muM,
                              slope = focal_mod$sims.list$betaM,
                              foc_season = "migration_low_evi")
# high evi
cond_mh<- getConditionalSurv(model = focal_mod, covar = eviq[2], 
                             intercept = focal_mod$sims.list$muM,
                             slope = focal_mod$sims.list$betaM,
                             foc_season = "migration_high_evi")

# Bind results
data_mig <- rbind(cond_ml, cond_mh) %>% 
  rename(name = season, 
         est = estimate,
         upper = ci_h,
         lower = ci_l) %>% 
  mutate(dur = (17/7)) 

mig_evi <- compareSurv(data=data_mig)

# Combine results and write out
evi_comb <- rbind(wint_evi, breed_evi, mig_evi)
write_csv(evi_comb, path = "output/comp_evi_risk.csv")
