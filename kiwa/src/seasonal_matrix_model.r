# Set up simple seasonal surviavl model (4 seasons)
library(popbio)
library(rjags)
library(jagsUI)
library(tidyverse)
library(ggthemes)

`%notin%` <- Negate(`%in%`)

# Define a cyclic transition matrix as a function
buildProjection <- function(Fb, Psb, Pbf, Pfw, Pws, ...){
  tmat <- matrix(
    c(
      0,   0,   0, Psb,
      Pbf+Fb,   0,   0,   0,
      0, Pfw,   0,   0,
      0,   0, Pws,   0
    ), nrow = 4, byrow = T)
  return(tmat)
}

# version with post-breeding
# buildProjection <- function(Fb, Psb, Pbf, Pfw, Pws, Pbp){
#   tmat <- matrix(
#     c(
#       0,   0,   0,   0,  Psb,
#       Pbp,  0,  0,   0,  0,
#       0,  Pbf+Fb,   0,   0,   0,
#       0,  0, Pfw,   0,   0,
#       0,  0,   0, Pws,   0
#     ), nrow = 5, byrow = T)
#   return(tmat)
# }

# start with Rockwell's parameters
Fb <- 2.75      # breeding season fecundity (pulled from Birds of the World...)
Pbf <- 0.963^4  # breeding to fall survival
Pfw <- 0.879  # fall to winter survival
Pws <- 0.977^6  # winter to spring survival
Psb <- 0.879  # spring to breeding survival
Pbp <- 1

#annual survival
Pbf*Pfw*Pws*Psb


tmat <- buildProjection(Fb = Fb, Psb = Psb, Pbf = Pbf, Pfw = Pfw, Pws = Pws, Pbp = Pbp)

# extract lambda
lambda(tmat)




# hard code a simulation
N1 <- c(100, 0, 0, 0)
(N2 <- tmat%*%N1)
(N3 <- tmat%*%N2)
(N4 <- tmat%*%N3)
(N5 <- tmat%*%N4)


# now update use our estimates where we can
Pbf <- 0.964^(130/7)
Pws <- 0.968^(189/7)
Psb <- 0.899^(17/7)

#annual survival
Pbf*Pfw*Pws*Psb


tmat2 <- buildProjection(Fb = Fb, Psb = Psb, Pbf = Pbf, Pfw = Pfw, Pws = Pws, Pbp = Pbp)
lambda(tmat2)

# hard code a simulation
N1 <- c(100, 0, 0, 0)
(N2 <- tmat2%*%N1)
(N3 <- tmat2%*%N2)
(N4 <- tmat2%*%N3)
(N5 <- tmat2%*%N4)

# Finally, we use the coeffs to estimate it...
load("output/model_w2017_noage_2.rdata")

# define a function to return prob. of survival based on an evi value
predictSurv <- function(focal_mod, evi, season, error = c(0.2, 0.8)){
  
  if(season %notin% c("migration", "breeding", "winter")){
    message("Season argument must be one of: migration, breeding, or winter...")
  }
  
  if(season == "migration"){
    migpost <- plogis(focal_mod$sims.list$muM+focal_mod$sims.list$betaM*evi)
    meanS <- mean(migpost)
    lciS <- quantile(migpost, probs = c(error[1]))
    hciS <- quantile(migpost, probs = c(error[2]))
  }
  
  if(season == "breeding"){
    breedpost <- plogis(focal_mod$sims.list$muB+focal_mod$sims.list$betaB*evi)
    meanS <- mean(breedpost)
    lciS <- quantile(breedpost, probs = c(error[1]))
    hciS <- quantile(breedpost, probs = c(error[2]))
  }
  
  if(season == "winter"){
    wintpost <- plogis(focal_mod$sims.list$muW+focal_mod$sims.list$betaW*evi)
    meanS <- mean(wintpost)
    lciS <- quantile(wintpost, probs = c(error[1]))
    hciS <- quantile(wintpost, probs = c(error[2]))
  }
  
  
  return(c("mean" = meanS, "lower" = lciS, "upper" = hciS))
}

# predictSurv(focal_mod = mod, evi = 0.5, season = "breeding")

#get evi vals to interrogate effect on lambda
load("output/kiwa_dat_JAGS.rdata")

ps <- seq(0, 1, by = 0.05)

# set up vector of observed quantiles of evi
evi_vec <- quantile(evi, probs = ps, na.rm = T)

# pull vectors of seasonal survivals across the vector evi vals
wint_surv_vec <- sapply(evi_vec, FUN = predictSurv, focal_mod = mod, season = "winter")
breed_surv_vec <- sapply(evi_vec, FUN = predictSurv, focal_mod = mod, season = "breeding")
mig_surv_vec <- sapply(evi_vec, FUN = predictSurv, focal_mod = mod, season = "migration")

# define function to get 
iterateLambdas <- function(wint, breed, mig){
  Fb <- 2.75       # breeding season fecundity (pulled from Birds of the World...)
  Pfw <- 0.879  # estimate from Rockwell et al. 2017
  Pbf <- breed^(130/7)
  Pws <- wint^(189/7)
  Psb <- mig^(17/7)
  Pbp <- 1
  tmat <- buildProjection(Fb = Fb, Psb = Psb, Pbf = Pbf, Pfw = Pfw, Pws = Pws, Pbp = Pbp)
  l <- lambda(tmat)
  return(l)
}

iterateLambdas(wint = wint_surv_vec[1], breed = breed_surv_vec[1], mig = mig_surv_vec[1])

mean_lambda <- mapply(FUN = iterateLambdas, wint_surv_vec[1,], breed_surv_vec[1,], mig_surv_vec[1,])
low_lambda <- mapply(FUN = iterateLambdas, wint_surv_vec[2,], breed_surv_vec[2,], mig_surv_vec[2,])
hi_lambda <- mapply(FUN = iterateLambdas, wint_surv_vec[3,], breed_surv_vec[3,], mig_surv_vec[3,])

plot_df <- tibble(mean = mean_lambda, 
                  low = low_lambda, 
                  high = hi_lambda, 
                  evi_quant = ps,
                  group = 1)


(lambda_plot <- ggplot(plot_df) +
    # geom_ribbon(aes(x = evi_quant, ymin = low, ymax = high), linewidth = 2, alpha = 0.2, lineend = "round")+
    geom_point(aes(x = evi_quant, y = mean), size = 2, se = F, color = "black")+
    geom_hline(aes(yintercept = 1), linetype = "dashed")+
    xlab("Quantile of Observed EVI")+
    ylab(expression(paste("Population Growth Rate (",lambda, ")")))+
    ylim(c(0.8, 1.05))+
    coord_fixed(10/2.5)+
    # ggtitle("KIWA population growth rate as a fucntion of winter EVI")+
    theme_few(base_size = 14)
)

ggsave(lambda_plot, file = "figures/fig4 - lambda_plot.pdf", width = 4, height = 4)
ggsave(lambda_plot, file = "figures/fig4 - lambda_plot.png", width = 4, height = 4)

(lambda_plot_alt <- ggplot(plot_df) +
    geom_point(aes(x = evi_quant, y = mean), size = 3)+
    geom_line(aes(x = evi_quant, y = mean, group = group), size = 1, linetype = 5)+
    # geom_errorbar(aes(x = evi_quant, ymin = low, ymax = high), width = 0.3, linewidth = 1)+
    geom_hline(aes(yintercept = 1), linetype = "dashed", alpha = 0.8)+
    xlab("Quantile of Observed EVI")+
    ylab(expression(paste("Population Growth Rate (",lambda, ")")))+
    ylim(c(0.85, 1.05))+
    # ggtitle("KIWA population growth rate as a fucntion of winter EVI")+
    theme_few(base_size = 14)
)

ggsave(lambda_plot_alt, file = "figures/fig4 - lambda_plot_alt.pdf", device = cairo_pdf(), width = 4, height = 3)
