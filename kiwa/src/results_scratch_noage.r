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


#extract model numbers

mods_nums <- c(1,2)
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

# Beta table

(betas <- data.frame(season = c("Winter", "Migration", "Breeding"),
                    mean = plogis(c(focal_mod$mean$betaW, focal_mod$mean$betaM, focal_mod$mean$betaB)),
                    cil = plogis(c(focal_mod$q2.5$betaW, focal_mod$q2.5$betaM, focal_mod$q2.5$betaB)),
                    cih = plogis(c(focal_mod$q97.5$betaW, focal_mod$q97.5$betaM, focal_mod$q97.5$betaB)),
                    p = c(sum(plogis(focal_mod$sims.list$betaW) > 0)/length(focal_mod$sims.list$betaW),
                          sum(plogis(focal_mod$sims.list$betaM) > 0)/length(focal_mod$sims.list$betaM),
                          sum(plogis(focal_mod$sims.list$betaB) > 0)/length(focal_mod$sims.list$betaB))))

evi_new <- seq(-1, 1, length.out = 100)

## Make Conditional Predictions ##

# Winter
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

# Migration
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

# Breeding
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

# Combine predictions
dat <- rbind(dat_W, dat_M, dat_B)

##-- Plotting  --##

# Define a color palette
pal  <- c("#547AA5", "#2C6E49", "#A44200")


#load data to extract covariate values lims
load("output/kiwa_dat_JAGS.rdata")

# Extract EVI
raws <- data.frame(evi = evi)
min_evi <- min(na.omit(evi))
max_evi <- max(na.omit(evi))
mean_evi <- mean(na.omit(evi))


# Marginal seasonal survival plot
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

## Get conditioonal pred at mean evi for each season ##

# Breeding
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

# Migration
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

# Winter
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

# Combine into single df
conds <- rbind(cond_b, cond_w, cond_m)


## Dot and whisker plot

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

# EXTRA ANALYSES NOT USED IN PAPER

# Get marginal slopes
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

ggsave(pheno_plot, filename = "output/pheno_plot.pdf", width = 8, height = 4)
ggsave(pheno_plot, filename = "output/pheno_plot.png", width = 8, height = 4)

