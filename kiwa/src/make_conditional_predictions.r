library(tidyverse)
library(glue)
library(jagsUI)
library(emmeans)
library(ggthemes)
# library(R2jags)
library(rjags)
library(tidybayes)
# options(max.print = 1000)

# get list of model rdata files
mods <- list.files("output", full.names = T, pattern = "noage")
# [18:19]

#extract model numbers
# mods_nums <- as.numeric(str_extract(mods, "[0-9]+"))
mods_nums <- c(1,2)
# stupid for loop I had to make b/c I did a bad job naming model objects upon creation...
for(i in 1:length(mods)){
  load(mods[i])
  assign(glue("mod{mods_nums[i]}"), mod)
}


#define focal mod for plotting
focal_mod <- mod2





#--  Get conditional posteriors --#

# Beta table

(betas <- data.frame(season = c("Winter", "Migration", "Breeding"),
                    mean = plogis(c(focal_mod$mean$betaW, focal_mod$mean$betaM, focal_mod$mean$betaB)),
                    cil = plogis(c(focal_mod$q2.5$betaW, focal_mod$q2.5$betaM, focal_mod$q2.5$betaB)),
                    cih = plogis(c(focal_mod$q97.5$betaW, focal_mod$q97.5$betaM, focal_mod$q97.5$betaB)),
                    p = c(sum(plogis(focal_mod$sims.list$betaW) > 0)/length(focal_mod$sims.list$betaW),
                          sum(plogis(focal_mod$sims.list$betaM) > 0)/length(focal_mod$sims.list$betaM),
                          sum(plogis(focal_mod$sims.list$betaB) > 0)/length(focal_mod$sims.list$betaB))))

evi_new <- seq(-1, 1, length.out = 100)

post <- list()
mean <- c()
quant_l <- c()
quant_h <- c()
tmp <- list()
age_class <- c("SY", "ASY")
for(i in 1:length(evi_new)){
  post[[i]] <- plogis(focal_mod$sims.list$muW+focal_mod$sims.list$betaW*evi_new[i])
  mean[i] <- mean(post[[i]])
  quant_l[i] <- quantile(post[[i]], probs = c(0.025))
  quant_h[i] <- quantile(post[[i]], probs = c(0.975))
  
}
dat_W <- data.frame(evi = evi_new,
                    mean = mean,
                    cil = quant_l,
                    cih = quant_h,
                    season = "Winter")


post <- list()
mean <- c()
quant_l <- c()
quant_h <- c()
tmp <- list()
for(i in 1:length(evi_new)){
  post[[i]] <- plogis(focal_mod$sims.list$muM+focal_mod$sims.list$betaM*evi_new[i])
  mean[i] <- mean(post[[i]])
  quant_l[i] <- quantile(post[[i]], probs = c(0.025))
  quant_h[i] <- quantile(post[[i]], probs = c(0.975))
  
}
dat_M <- data.frame(evi = evi_new,
                    mean = mean,
                    cil = quant_l,
                    cih = quant_h, 
                    season = "Migration")


post <- list()
mean <- c()
quant_l <- c()
quant_h <- c()
tmp <- list()
for(i in 1:length(evi_new)){
  post[[i]] <- plogis(focal_mod$sims.list$muB+focal_mod$sims.list$betaB*evi_new[i])
  mean[i] <- mean(post[[i]])
  quant_l[i] <- quantile(post[[i]], probs = c(0.025))
  quant_h[i] <- quantile(post[[i]], probs = c(0.975))
  # 
}
dat_B <- data.frame(evi = evi_new,
                    mean = mean,
                    cil = quant_l,
                    cih = quant_h, 
                    season = "Breeding")

dat <- rbind(dat_W, dat_M, dat_B)
write_rds(dat, "figures/kiwa_data/evi_pred.rds")

#load data to extratc covar lims
load("output/kiwa_dat_JAGS.rdata")
raws <- data.frame(evi = evi)
write_rds(raws, "figures/kiwa_data/evi_dat.rds")

