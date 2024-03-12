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


# Extract DIC
DICs <- c()
for(i in 1:2){
  m <- get(glue("mod{mods_nums[i]}"))
  DICs[i] <- m$DIC
}

dic_df <- data.frame(DIC = DICs,
                     model_no = c(1,2))

# Add DIC to model descriptions file
mod_des <- read.csv("ctfs/model_control.csv") %>% 
  left_join(dic_df)


# Print (sorted by DIC)
mod_des %>% 
  arrange(DIC)


#define focal mod for plotting
focal_mod <- mod2

# Print model summary
options(max.print = 10000)
focal_mod
options(max.print = 1000)
summary(focal_mod)

jagsUI::traceplot(focal_mod, parameters = c(
  "muW", "betaW",
  "muM", "betaM",
  "muB", "betaB"
))

focal_mod %>%
  spread_draws(muW, muM, muB) %>%
  ggplot(aes(y = fct_rev(condition), x = condition_mean)) +
  stat_halfeye(.width = c(.90, .5))

densityplot(focal_mod, parameters = c("muM", "betaM", "muW", "betaW", "mean.phiB"))

#----------------------------------#
# get prop posterior > 0 for slope params




#----------------------------------#

# Beta table

(betas <- data.frame(season = c("Winter", "Migration", "Breeding"),
                    mean = plogis(c(focal_mod$mean$betaW, focal_mod$mean$betaM, focal_mod$mean$betaB)),
                    cil = plogis(c(focal_mod$q2.5$betaW, focal_mod$q2.5$betaM, focal_mod$q2.5$betaB)),
                    cih = plogis(c(focal_mod$q97.5$betaW, focal_mod$q97.5$betaM, focal_mod$q97.5$betaB)),
                    p = c(sum(plogis(focal_mod$sims.list$betaW) > 0)/length(focal_mod$sims.list$betaW),
                          sum(plogis(focal_mod$sims.list$betaM) > 0)/length(focal_mod$sims.list$betaM),
                          sum(plogis(focal_mod$sims.list$betaB) > 0)/length(focal_mod$sims.list$betaB))))

evi_new <- seq(-1, 1, length.out = 100)

# TODO: make this a function

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

# dat_B <- data.frame(evi = rep(evi_new, 2),
#                     mean = rep(focal_mod$mean$mean.phiB, 100),
#                     cil = rep(focal_mod$q2.5$mean.phiB, 100),
#                     cih = rep(focal_mod$q97.5$mean.phiB, 100),
#                     season = rep("Breeding", 100),
#                     age = rep(age_class, 100))
# 

dat <- rbind(dat_W, dat_M, dat_B)
pal  <- c("#547AA5", "#2C6E49", "#A44200")



#load data to extratc covar lims
load("output/kiwa_dat_JAGS.rdata")
raws <- data.frame(evi = evi)

min_evi <- min(na.omit(evi))
max_evi <- max(na.omit(evi))
mean_evi <- mean(na.omit(evi))



(marg <- ggplot()+
    geom_line(data = dat, aes(x=evi, y = mean, color = season), size = 1.3) +
    geom_ribbon(data = dat, aes(x=evi, ymin = cil, ymax = cih, fill = season), alpha = 0.2) +
    scale_color_manual(values = pal, name = "Season") +
    scale_fill_manual(values = pal, name = "Season") +
    geom_rug(data = raws, aes(x = evi), sides = "b") +
    ylab("Weekly Survival") +
    xlab("EVI") +
    theme_minimal()+
    xlim(c(0,max_evi+0.1))+
    ylim(c(0.75,1))+
    # ylim(c(0.7,1))+
    # facet_wrap(~age)+
    theme_classic(base_size = 18) +
    theme(panel.spacing = unit(2, "lines")))

ggsave(marg, filename = "output/survival_plot.pdf")
ggsave(marg, filename = "output/survival_plot.png", width = 8,
       height = 5)

#get conditioonal pred at mean evi
# breeding
m <- c()
quant_l <- c()
quant_h <- c()
# age_l <- c()
post <- plogis(focal_mod$sims.list$muB+focal_mod$sims.list$betaB*mean_evi)
m <- median(post)
quant_l <- quantile(post, probs = c(0.025))
quant_h <- quantile(post, probs = c(0.975))
# age_l <- age_class
season <- "breeding"

cond_b  <- data.frame(estimate = m,
                      ci_l = quant_l,
                      ci_h = quant_h,
                      # age = age_l,
                      season = season)

m <- c()
quant_l <- c()
quant_h <- c()
age_l <- c()
  post <- plogis(focal_mod$sims.list$muW+focal_mod$sims.list$betaW*mean_evi)
  m <- median(post)
  quant_l <- quantile(post, probs = c(0.025))
  quant_h <- quantile(post, probs = c(0.975))
  # age_l <- age_class
  season <- "winter"

cond_w  <- data.frame(estimate = m,
                      ci_l = quant_l,
                      ci_h = quant_h,
                      season = season)

m <- c()
quant_l <- c()
quant_h <- c()
  post <- plogis(focal_mod$sims.list$muM+focal_mod$sims.list$betaM*mean_evi)
  m <- median(post)
  quant_l <- quantile(post, probs = c(0.025))
  quant_h <- quantile(post, probs = c(0.975))
  # age_l <- age_class
  season <- "migration"


cond_m  <- data.frame(estimate = m,
                      ci_l = quant_l,
                      ci_h = quant_h,
                      # age = age_l,
                      season = season)

conds <- rbind(cond_b, cond_w, cond_m)

write_csv(conds, "output/files_for_nate/conditional_effects_noage.csv")

pal  <- c("#547AA5", "#2C6E49", "#A44200")

(marg <- ggplot()+
    geom_line(data = dat, aes(x=evi, y = mean, color = season)) +
    geom_ribbon(data = dat, aes(x=evi, ymin = cil, ymax = cih, fill = season), alpha = 0.1) +
    scale_color_manual(values = pal, name = "Season") +
    scale_fill_manual(values = pal, name = "Season") +
    geom_rug(data = raws, aes(x = evi), sides = "b") +
    
    xlab("EVI") +
    theme_minimal()+
    xlim(c(0,max_evi+0.1))+
    # ylim(c(0.7,1))+
    # facet_wrap(~age)+
    theme(panel.spacing = unit(2, "lines")))

ggsave(marg, filename = "output/files_for_nate/survival_plot.pdf")
ggsave(marg, filename = "output/files_for_nate/survival_plot.png", width = 8,
       height = 5)

# save top model
save(focal_mod, file = "output/files_for_nate/top_mod.rdata")


# Dot and whisker plot

dot_df <- data_frame(season = c("Winter", "Migration", "Breeding"),
                     # age = c(rep(age_class, 3)),
                     mu = c(mod1$mean$mean.phiW, 
                            mod1$mean$mean.phiM, 
                            mod1$mean$mean.phiB),
                     cil = c(mod1$q2.5$mean.phiW,
                             mod1$q2.5$mean.phiM,
                             mod1$q2.5$mean.phiB),
                     cih = c(mod1$q97.5$mean.phiW,
                             mod1$q97.5$mean.phiM,
                             mod1$q97.5$mean.phiB))


(dot_marg <- ggplot(dot_df) +
    geom_pointrange(aes(x = season, y=mu, ymin = cil, ymax=cih,
                        color = season),
                    size = 1, linewidth = 1.2,
                    position=position_dodge(width=1)) +
    scale_color_manual(values = pal, name = "Season") +
    ylab("Weekly Survival") +
    theme_minimal()
)

slopes_df <- data_frame(season = c("Winter", "Migration", "Breeding"),
                        # age = c(rep(age_class, 3)),
                        mu = c(plogis(mod15$mean$betaW), 
                               plogis(mod15$mean$betaM), 
                               plogis(mod15$mean$betaB)),
                        cil = c(plogis(mod15$q2.5$betaW),
                                plogis(mod15$q2.5$betaM),
                                plogis(mod15$q2.5$betaB)),
                        cih = c(plogis(mod15$q97.5$betaW),
                                plogis(mod15$q97.5$betaM),
                                plogis(mod15$q97.5$betaB)))
(int_marg <- ggplot(dot_df) +
    geom_pointrange(aes(x = season, group = age, y=mu, ymin = cil, ymax=cih,
                        color = season, shape = age),
                    position=position_dodge(width=1)) +
    theme_minimal()
)



# Transition curves
pheno_w2m <- data.frame(mean = focal_mod$mean$depWM,
                        cil = focal_mod$q2.5$depWM,
                        cih = focal_mod$q97.5$depWM,
                        trans = "Winter to Migration", 
                        week = 1:17)

pheno_m2b <- data.frame(mean = focal_mod$mean$depMB,
                        cil = focal_mod$q2.5$depMB,
                        cih = focal_mod$q97.5$depMB,
                        trans = "Migration to Breeding",
                        week = 1:17)

pheno <- rbind(pheno_w2m, pheno_m2b)

pheno_plot <- ggplot(pheno)+
  geom_point(aes(x=week, y = mean, color = trans), 
             position = position_dodge(0.5), size = 1)+
  geom_errorbar(aes(x=week, ymin = cil, ymax = cih, color = trans), 
                position = position_dodge(0.5), size = 0.6, width = 0.8) +
  ylab("Transition Probability")+
  xlab("Study Week") +
  scale_color_manual(values = pal[1:3], name = "Transition")+
  theme_minimal()+
  NULL
  pheno_plot

ggsave(pheno_plot, filename = "output/files_for_nate/pheno_plot.pdf", width = 8, height = 4)
ggsave(pheno_plot, filename = "output/files_for_nate/pheno_plot.png", width = 8, height = 4)

zmat <- focal_mod$q50$z


# Hyperparams for the Weibull params
# gWM ~ dunif(0, 17)
gMB <- 8.8
# gMB ~ dunif(0, 17)

# kWM ~ dgamma(1, 1)
kMB <- 3.2
# kMB ~ dgamma(1, 1)


# check priors on transition
depMB <- c()
for(t in 6:17) {
  depMB[t] <- 1-exp(-(t/gMB)^kMB + ((t-1)/gMB)^kMB)
}

ggplot(pheno)+
  geom_point(aes(x=week, y = mean, color = trans))+
  geom_point(data = data.frame(prior = depMB, x = 1:17), aes(x=x, y = prior))+
  geom_errorbar(aes(x=week, ymin = cil, ymax = cih, color = trans), alpha = 0.25) +
  ylab("Cummulative Transition Probability")+
  theme_minimal()


#----------------------------------#
# SCRATCH/DEPRECATED #


load("output/kiwa_dat_JAGS.rdata")

ggplot(data.frame(evi = evi))+
  geom_density(aes(x=evi))

min(na.omit(evi))
max(na.omit(evi))
sum(na.omit(evi) > 0)


muW_SY <- focal_mod$mean$muW[1]
betaW_SY <- focal_mod$mean$betaW[1]
muW_SY_ucl <- focal_mod$q2.5$muW[1]
betaW_SY_ucl <- focal_mod$q2.5$betaW[1]
muW_SY_lcl <- focal_mod$q97.5$muW[1]
betaW_SY_lcl <- focal_mod$q97.5$betaW[1]

muW_AHY <- focal_mod$mean$muW[2]
betaW_AHY <- focal_mod$mean$betaW[2]
muW_AHY_ucl <- focal_mod$q2.5$muW[2]
betaW_AHY_ucl <- focal_mod$q2.5$betaW[2]
muW_AHY_lcl <- focal_mod$q97.5$muW[2]
betaW_AHY_lcl <- focal_mod$q97.5$betaW[2]

muM_SY <- focal_mod$mean$muM[1]
betaM_SY <- focal_mod$mean$betaM[1]
muM_SY_ucl <- focal_mod$q2.5$muM[1]
betaM_SY_ucl <- focal_mod$q2.5$betaM[1]
muM_SY_lcl <- focal_mod$q97.5$muM[1]
betaM_SY_lcl <- focal_mod$q97.5$betaM[1]

muM_AHY <- focal_mod$mean$muM[2]
betaM_AHY <- focal_mod$mean$betaM[2]
muM_AHY_ucl <- focal_mod$q2.5$muM[2]
betaM_AHY_ucl <- focal_mod$q2.5$betaM[2]
muM_AHY_lcl <- focal_mod$q97.5$muM[2]
betaM_AHY_lcl <- focal_mod$q97.5$betaM[2]

phiB_SY <- focal_mod$mean$mean.phiB[1]
phiB_SY_ucl <- focal_mod$q2.5$mean.phiB[1]
phiB_SY_lcl <- focal_mod$q2.5$mean.phiB[1]

phiB_AHY <- focal_mod$mean$mean.phiB[2]
phiB_AHY_ucl <- focal_mod$q2.5$mean.phiB[2]
phiB_AHY_lcl <- focal_mod$q2.5$mean.phiB[2]


evi <- seq(-1, 1, length.out = 100)

W_pred_SY <- data.frame(
  evi = evi,
  mean = plogis(muW_SY+betaW_SY*evi),
  ucl = plogis(muW_SY_ucl+betaW_SY_ucl*evi),
  lcl = plogis(muW_SY_lcl+betaW_SY_lcl*evi),
  season = rep("winter", 100),
  age = rep("SY", 100)
)

W_pred_AHY <- data.frame(
  evi = evi,
  mean = plogis(muW_AHY+betaW_AHY*evi),
  ucl = plogis(muW_AHY_ucl+betaW_AHY_ucl*evi),
  lcl = plogis(muW_AHY_lcl+betaW_AHY_lcl*evi),
  season = rep("winter", 100),
  age = rep("AHY", 100)
)

M_pred_SY <- data.frame(
  evi = evi,
  mean = plogis(muM_SY+betaM_SY*evi),
  ucl = plogis(muM_SY_ucl+betaM_SY_ucl*evi),
  lcl = plogis(muM_SY_lcl+betaM_SY_lcl*evi),
  season = rep("migration", 100),
  age = rep("SY", 100)
)

M_pred_AHY <- data.frame(
  evi = evi,
  mean = plogis(muM_AHY+betaM_AHY*evi),
  ucl = plogis(muM_AHY_ucl+betaM_AHY_ucl*evi),
  lcl = plogis(muM_AHY_lcl+betaM_AHY_lcl*evi),
  season = rep("migration", 100),
  age = rep("AHY", 100)
)

B_pred_SY <- data.frame(
  evi = evi,
  mean = phiB_SY,
  ucl = phiB_SY_ucl,
  lcl = phiB_SY_lcl,
  season = rep("breeding", 100),
  age = rep("SY", 100)
)

B_pred_AHY <- data.frame(
  evi = evi,
  mean = phiB_AHY,
  ucl = phiB_AHY_ucl,
  lcl = phiB_AHY_lcl,
  season = rep("breeding", 100),
  age = rep("AHY", 100)
)

pred_df <- rbind(W_pred_SY, W_pred_AHY, M_pred_SY, M_pred_AHY, B_pred_SY, B_pred_AHY)


ggplot(pred_df) +
  geom_line(aes(x=evi, y = mean, group = season, color = season))+
  # geom_ribbon(aes(x = evi, ymin=lcl, ymax=ucl, fill = season, group = season), alpha = 0.25) +
  ylim(0,1)+
  facet_wrap(~age) +
  theme_minimal()

migs <- c()
for(i in 1:17) {
  v <- zmat[,i]
  twos <- sum(v == 2)
  fives <- sum(v == 5)
  tot <- twos + fives
  migs[i] <- tot 
}
migs

breeds <- c()
for(i in 1:17) {
  v <- zmat[,i]
  twos <- sum(v == 3)
  fives <- sum(v == 6)
  tot <- twos + fives
  breeds[i] <- tot 
}
breeds

