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


btbw <- data.frame(name = c("spring", "fall", "summer", "winter"),
                     est = c(0.967, 0.970, 0.994, 0.997),
                     upper = c(0.998, 0.999, 0.999, 0.992),
                     lower = c(0.932, 0.930, 0.987, 0.999),
                     dur = c(1,1,1,1))

# Compare one weeks
sp_fa_df <- btbw[c(1,2),]
sp_fa_out <- compareSurv(sp_fa_df)

su_wi_df <- btbw[c(3,4),]
su_wi_out <- compareSurv(su_wi_df)

sp_wi_df <- btbw[c(1,4),]
sp_wi_out <- compareSurv(sp_wi_df)

fa_wi_df <- btbw[c(2,4),]
fa_wi_out <- compareSurv(fa_wi_df)

sp_su_df <- btbw[c(1,3),]
sp_su_out <- compareSurv(sp_su_df)

fa_su_df <- btbw[c(2,4),]
fa_su_out <- compareSurv(fa_su_df)

# 