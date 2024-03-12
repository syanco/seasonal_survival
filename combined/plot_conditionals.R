#################################################
####                                         ####
####   Plot predicted conditional effects.   ####
####                                         ####
#################################################

####---- Inits ----####
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggthemes)

# pal  <- c("#547AA5", "#2C6E49", "#A44200")
# pal <- c("#003366", "#99CC00", "#FF3300")#opt1
# pal <- c("#4169E1", "#50C878", "#FF6347")# opt 2
# pal <- c("#4575B4", "#99CC00", "#D73027")# opt 3
pal <- c("#003366", "#4CAF50", "#FF3300") #opt 4


####---- Load Data ----####

#- BTBW -#
cov_df <- readRDS("figures/btbw_data/cov_df.rds")
pred_soi <- readRDS("figures/btbw_data/pred_soi.rds")%>% 
  mutate(season = case_when(Season == "Migration" ~ "Spring migration",
                            Season == "Winter" ~ "Non-breeding",
                            T ~ Season),
         season = fct_relevel(season, "Non-breeding", "Spring migration"))
pred_bseas <- readRDS("figures/btbw_data/pred_bseas.rds")%>% 
  mutate(season = case_when(Season == "Autumn migration" ~ "Fall migration",
                            Season == "Winter" ~ "Non-breeding",
                            T ~ Season),
         season = fct_relevel(season, "Breeding", "Fall migration"))

#- KIWA -#
evi_df <- read_rds("figures/kiwa_data/evi_dat.rds")
pred_evi <- read_rds("figures/kiwa_data/evi_pred.rds") %>% 
  mutate(season = case_when(season == "Migration" ~ "Spring migration",
                            season == "Winter" ~ "Non-breeding",
                            T ~ season),
         season = fct_relevel(season, "Non-breeding", "Spring migration"))

arrow_data <- tibble(x = 0.25,
                     xend = 0.75,
                     y = 0.76, 
                     yend = 0.76)
# Sample data for horizontal arrow and labels
arrow_data <- data.frame(
  x = c(0.3, 0.7),   # X-coordinates for arrow start points
  xend = c(0.7, 0.3), # X-coordinates for arrow end points
  y = c(0.77, 0.77),       # Y-coordinates for arrow start points
  yend = c(0.77, 0.77),     # Y-coordinates for arrow end points
  lab_x = c(0.2, 0.8),
  label = c("Dry", "Wet")  # Labels for each segment
)

####---- Process Data ----###

# get covar bounds for plotting

# EVI
min_evi <- min(na.omit(evi_df$evi))
max_evi <- max(na.omit(evi_df$evi))
mean_evi <- mean(na.omit(evi_df$evi))

# SOI
min_soi <- min(na.omit(cov_df$SOI))
max_soi <- max(na.omit(cov_df$SOI))
mean_Soi <- mean(na.omit(cov_df$SOI))


####---- Make Figs ----####


##-- Fig 2 --##

(kiwa_cond <- ggplot()+
   geom_line(data = pred_evi, aes(x=evi, y = mean, color = season), size = 1, alpha = 1) +
   geom_ribbon(data = pred_evi, aes(x=evi, ymin = cil, ymax = cih, fill = season), alpha = 0.1) +
   # scale_linetype_manual(values = c(2,3,4), name = "Season")+ 
   scale_color_manual(values = pal, name = "Season") +
   scale_fill_manual(values = pal, name = "Season") +
   geom_rug(data = evi_df, aes(x = evi), sides = "b") +
   ylab("Weekly Survival") +
   xlab("Winter Enhanced Vegetation Index") +
   ggtitle("Kirtland's Warbler") +
   # theme_minimal()+
   xlim(c(0,max_evi+0.1))+
   ylim(c(0.75,1))+
   # ylim(c(0.7,1))+
   # facet_wrap(~age)+
   theme_few(base_size = 14) +
   # theme(panel.spacing = unit(2, "lines"),
   #       plot.title = element_text(hjust = 0.5))+
   coord_fixed(ratio = 4)+
   # geom_segment(data = arrow_data, aes(x = x, xend = xend, y = y, yend = yend),
   #              arrow = arrow(type = "closed", angle = 20), size = 1) +
   # geom_text(data = arrow_data, aes(x = lab_x, y = y, label = label),
   #           vjust = 0, size = 4)+
   NULL
 )


(btbw_cond <- ggplot()+
    geom_line(data = pred_soi, aes(x=SOI, y = Mean, color = season), size = 1) +
    geom_ribbon(data = pred_soi, aes(x=SOI, ymin = LCI.95, ymax = UCI.95, fill = season), alpha = 0.1) +
    scale_linetype_manual(values = c(1,2,4), name = "Season")+ 
    scale_color_manual(values = pal, name = "Season") +
    scale_fill_manual(values = pal, name = "Season") +
    geom_rug(data = cov_df, aes(x = SOI), sides = "b") +
    ylab("") +
    xlab("Southern Oscillation Index") +
    ggtitle("Black-throated Blue Warbler") +
    # theme_minimal()+
    xlim(c(min_soi-0.1,max_soi+0.1))+
    ylim(c(0.85,1))+
    # ylim(c(0.7,1))+
    # facet_wrap(~age)+
    theme_few(base_size = 14) +
    theme(panel.spacing = unit(.5, "lines"),
          plot.title = element_text(hjust = 0.5))+
    facet_wrap(~Sex, nrow = 2, strip.position = "right")+
    coord_fixed(ratio = 7)+
    NULL
  )

# (kiwa_cond/btbw_cond) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
(comb_plot <- (kiwa_cond|btbw_cond) + plot_layout(guides = "collect") & theme(legend.position = "bottom") & 
    theme(plot.margin = margin(10, 10, 0, 1)))

ggsave(comb_plot, filename = "figures/fig2 - survival_plot.pdf", width = 8, height = 5)
ggsave(comb_plot, filename = "figures/fig2 - survival_plot.png", width = 8, height = 5)

##-- Fig 3 --##

(btbw_bseas_cond <- ggplot()+
    geom_line(data = pred_bseas, aes(x=BSEAS, y = Mean, color = season), size = 1.3) +
    geom_ribbon(data = pred_bseas, aes(x=BSEAS, ymin = LCI.95, ymax = UCI.95, fill = season), alpha = 0.2) +
    scale_linetype_manual(values = c(1,2,4), name = "Season")+ 
    scale_color_manual(values = pal, name = "Season") +
    scale_fill_manual(values = pal, name = "Season") +
    geom_rug(data = cov_df, aes(x = bseas), sides = "b") +
    ylab("Weekly Survival") +
    xlab("Standardized Breeding Season Length") +
    # theme_minimal()+
    xlim(c(-1, 1))+
    ylim(c(0.85,1))+
    # ylim(c(0.7,1))+
    # facet_wrap(~age)+
    theme_few(base_size = 14) +
    theme(panel.spacing = unit(.5, "lines"),
          legend.position = "bottom")+
    facet_wrap(~Sex, nrow = 1))

ggsave(btbw_bseas_cond, filename = "figures/fig3 - bseas_survival_plot.pdf", width = 6, height = 4)
ggsave(btbw_bseas_cond, filename = "figures/fig3 - bseas_survival_plot.png", width = 6, height = 4)


####---- Effects Ribbon ----####


kw_w_e <- 0.651
kw_s_e <- 0.620
kw_b_e <- 0.620

kw_w_d <- 189
kw_s_d <- 17
kw_b_d <- 130
tot <- 189+17+130
kw_dat <- data_frame(day = 1:(kw_w_d+kw_s_d+kw_b_d),
                     season = case_when(day %in% 1:kw_w_d ~ "Non-Breeding",
                                        day %in% (kw_w_d+1):(kw_w_d+kw_s_d) ~ "Migration",
                                        day %in% (kw_w_d+kw_s_d+1):(kw_w_d+kw_s_d+kw_b_d) ~ "Breeding"),
                     eff = c(rep(kw_w_e, kw_w_d), rep(kw_s_e, kw_s_d), rep(kw_b_e, kw_b_d)))

ggplot(kw_dat) +
  geom_line(aes(x = day, y = eff))+
  scale_color_manual(name = "Season", values = pal)+
  geom_hline(aes(yintercept = 0))+
  theme_void()+
  ylim(c(-1,1))+
  NULL


bt_w_e <- 0.0
bt_s_e <- 0.92
bt_b_e <- -1.47

bt_w_cil <- -2.88
bt_s_cil <- 0.02
bt_b_cil <- -3.34

bt_w_cih <- 1.24
bt_s_cih <- 2.04
bt_b_cih <- 0.34

bt_w_d <- 189
bt_s_d <- 17
bt_b_d <- 130

bt_dat <- data_frame(day = 1:(bt_w_d+bt_s_d+bt_b_d),
                     season = case_when(day %in% 1:bt_w_d ~ "Non-Breeding",
                                        day %in% (bt_w_d+1):(bt_w_d+bt_s_d) ~ "Migration",
                                        day %in% (bt_w_d+bt_s_d+1):(bt_w_d+bt_s_d+bt_b_d) ~ "Breeding"),
                     eff = c(rep(bt_w_e, bt_w_d), rep(bt_s_e, bt_s_d), rep(bt_b_e, bt_b_d)),
                     cil = c(rep(bt_w_cil, bt_w_d), rep(bt_s_cil, bt_s_d), rep(bt_b_cil, bt_b_d)),
                     cih = c(rep(bt_w_cih, bt_w_d), rep(bt_s_cih, bt_s_d), rep(bt_b_cih, bt_b_d)),
                     plt_min = case_when(eff >=0 ~ 0.1,
                                         T ~ eff),
                     plt_max = case_when(eff < 0 ~ 0, 
                                         T ~ eff))

ggplot(bt_dat) +
  # geom_ribbon(aes(x=day, ymin = cil, ymax = cih), alpha = 0.4)+
  # geom_line(aes(x = day, y = eff, color = season, group = 1), linewidth = 3)+
  geom_ribbon(aes(x = day, ymin = plt_min, ymax = plt_max), fill = "darkorange")+
  geom_hline(aes(yintercept = 0), linetype ="dashed")+
  # geom_ribbon(aes(x=day,ymin = -2, ymax = 2, fill = season), alpha = 0.2)+
  scale_fill_manual(name = "Season", values = pal)+
  theme_void()+
  ylim(c(-2,2))+
  NULL

