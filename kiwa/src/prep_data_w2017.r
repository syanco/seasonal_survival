####################################
####                            ####
####        Prep KIWA Data      ####
####      Scott Yanco, PhD      ####
####    scott.yanco@yale.edu    ####
####                            ####
####################################


####---- Description ----####

# This script executes preps all data needed to fit robust design multi-state 
# Cormack-Jolly-Seber (CJS) survival models for Kirtland's Warbler (KIWA).  
# Process include cleaning, and reshaping capture history data, making derived 
# capture history objects (e.g., for primary and secondary capture occasiosn),
# as well as covariates.  The script writes these and saves them to a .Rdata 
# file that contains all objects needed to fit any of the proposed models in
# subsequent scripts.


####---- Initialization ----####

# Libraries
library(tidyverse)
library(lubridate)
library(glue)
library(assertthat)

# External Source Files
source("src/kiwa_funs.r")


####---- Load Data ----####

# MOTUS Tracking dataset
load("data/KW_Survival_With_Age July 27 2021.rda")

# Environmental annotations
annos <- read_csv("output/kiwa_annos.csv") %>% 
  mutate(evi = case_when(`value_landsat8-evi-250-30` < -1 ~ NA_real_,
                         `value_landsat8-evi-250-30` > 1 ~ NA_real_,
                         TRUE ~ `value_landsat8-evi-250-30`))


####---- Process Data ----####


  #####-- Create Multi-state Cap History --#####

# Assign states based on detection tower (e.g., MI towers = breeding)
kiwa_dat <- KW_Summary %>% 
  group_by(motusTagID) %>% # group by tag ID
  filter(!duplicated(ts_Yday)) %>%  # remove duplicates within a single day
  ungroup()%>% # remove grouping
  mutate(year = year(ts)) %>% #create year var
  # filter(year != 2017) %>% 
  mutate(ch = 1) %>% #create single state cap history var [depricated]
  mutate(ch_ms = case_when( #create seasonal state cap history
    Period == "Tagged" ~ 1,
    Period == "Winter" ~ 1,
    Period == "Migration" ~ 2,
    Period == "Breeding" ~ 3,
    Period == "Died" ~ 4,
    #observations after breeding will get captured as breeding state for the 
    #purposes of feeding the model known states - these obs get censored below
    Period == "Breeding_Departure" ~ 3, 
    Period == "Fall_Migration" ~ 3))

# Convert to wide format
kiwa_wide_ch <- kiwa_dat %>% 
  ungroup() %>% #remove any remaining grouping
  #get full combination of tag Ids and sample occasions
  tidyr::expand(motusTagID, ts_Yday = full_seq(kiwa_dat$ts_Yday, 1)) %>% 
  left_join(kiwa_dat) %>% #join thisback up with the full dataset to retain vars 
  mutate(ch_ms = replace_na(ch_ms, 4)) %>% #replace  NAs with 4 (=not seen)
  #pivot to wide format for CJS model
  pivot_wider(id_cols = motusTagID, names_from = ts_Yday, values_from = ch_ms, 
              names_prefix = "O")
  
#hardcode remove erroneous ind
kiwa_wide_ch <- kiwa_wide_ch[kiwa_wide_ch$motusTagID != 29302,]



  #####-- Extract EVI Covariate --#####
kiwa_evi_key <- annos %>% 
  filter(!duplicated(motusTagID)) %>% # remove duplicates (shouldn't be any...)
  select(motusTagID, `value_landsat8-evi-250-30`) %>% #extract ID and variable
  right_join(kiwa_wide_ch) # join to CH data to match order

# Check that order matches the CH
assert_that(sum(kiwa_evi_key$motusTagID != kiwa_wide_ch$motusTagID) == 0)

# Extract as numeric object for JAGS
evi <- kiwa_evi_key$`value_landsat8-evi-250-30`


  #####-- Prepare CH data for JAGS --#####

#extract encounter history matrix
encounter0 <- as.matrix(kiwa_wide_ch)
encounter0 <- encounter0[,2:dim(encounter0)[2]]

#right censor for end of breeding season and exclude ID col
encounter <- encounter0[,1:119]

# Declare constants
# NOTE: hardcoded to this specific dataset!!!
n.primary <- 17      # number of primary occasions
n.secondary <- rep(7, 17)     # vector of secondary occasions
n.ind <- nrow(encounter) # number of birds
# Create index of secondary occasions
index <- list(1:7, # occasion index
              8:14,
              15:21,
              22:28,
              29:35,
              36:42,
              43:49,
              50:56,
              57:63,
              64:70,
              71:77,
              78:84,
              85:91,
              92:98,
              99:105,
              106:112,
              113:119)

# Calc number of individuals caught per primary period
caught <- rep(NA, n.primary) # init vector
for (i in 1:n.primary){
  tmp <- encounter[,index[[i]]]
  #since 28s would be all "4"s (i.e. no detection, count rows (inds) != 28) 
  caught[i] <- nrow(tmp[rowSums(tmp)!=28,])
}

# Make 3d array of encounter (inds X primaries X secondaries)
obs <- array(NA, dim = c(n.ind, n.primary, max(n.secondary))) # init array
for (i in 1:n.primary){
  obs[,i,1:n.secondary[i]] <- encounter[,index[[i]]]
}

# Check that dimensions are as expected
# when 2017 not included change 136 to 78
assert_that(sum(dim(obs) == c(136, 17, 7))==3) # check dims

# Make cap history matrix for primary periods only
ch <- matrix(NA, n.ind, n.primary) # init matrix
for (i in 1:n.ind){
  for (t in 1:n.primary){
    # print(any(obs[i,t,1:n.secondary[t]] != 4))
    ifelse(any(obs[i,t,1:n.secondary[t]] != 4), # if observed at all during primary
           # grab the highest state observed (excluding 4s for unobs)
           ch[i,t] <- max(obs[i,t,1:n.secondary[t]][obs[i,t,1:n.secondary[t]] != 4]), 
           ch[i,t] <- 4) # otherwise record a 4 (not seen)
  }
}

# Make modified secondary data and control sequences as in Rieke et al 
# (but do so by season)

#- WINTER

testW <- matrix(NA, n.ind, n.primary)
for (i in 1:nrow(testW)){
  for (j in 1:ncol(testW)){
    if(2 %notin% obs[i,j,] & 3 %notin% obs[i,j,]) {
      testW[i,j] <- sum(obs[i,j,]==1, na.rm = TRUE)
    }else{
      testW[i,j] <- 0
    }
  }
}

seenW <- array(NA, c(n.ind, n.primary, max(n.secondary)))
missedW <- array(NA, c(n.ind, n.primary, max(n.secondary)))

for (i in 1:nrow(testW)){
  for (t in 1:ncol(testW)){
    for (j in 1:n.secondary[t]){
      if(testW[i,t] > 1 & obs[i,t,j] == 1){seenW[i,t,j] <- 1}
      if(testW[i,t] >= 1 & obs[i,t,j] == 4){missedW[i,t,j] <- 1}
    }
  }
}

yesW <- matrix(NA, n.primary, max(n.secondary))
noW <- matrix(NA, n.primary, max(n.secondary))

for (i in 1:nrow(yesW)){
  for (j in 1:ncol(yesW)){
    yesW[i,j] <- sum(seenW[,i,j], na.rm = TRUE)
    noW[i,j] <- sum(missedW[,i,j], na.rm = TRUE)
  }
}

totalW <- yesW + noW
totalW

# MIGRATION
testM <- matrix(NA, n.ind, n.primary)
for (i in 1:nrow(testM)){
  for (j in 1:ncol(testM)){
    if(2 %in% obs[i,j,] & 3 %notin% obs[i,j,]) {
      testM[i,j] <- sum(obs[i,j,]==2 | obs[i,j,]==1, na.rm = TRUE)
    }else{
      testM[i,j] <- 0
    }
  }
}

seenM <- array(NA, c(n.ind, n.primary, max(n.secondary)))
missedM <- array(NA, c(n.ind, n.primary, max(n.secondary)))

for (i in 1:nrow(testM)){
  for (t in 1:ncol(testM)){
    for (j in 1:n.secondary[t]){
      if(testM[i,t] > 1 & obs[i,t,j] <= 2){seenM[i,t,j] <- 1} #1 or 2?
      if(testM[i,t] >= 1 & obs[i,t,j] == 4){missedM[i,t,j] <- 1}
    }
  }
}

yesM <- matrix(NA, n.primary, max(n.secondary))
noM <- matrix(NA, n.primary, max(n.secondary))

for (i in 1:nrow(yesM)){
  for (j in 1:ncol(yesM)){
    yesM[i,j] <- sum(seenM[,i,j], na.rm = TRUE)
    noM[i,j] <- sum(missedM[,i,j], na.rm = TRUE)
  }
}

totalM <- yesM + noM

# BREEDING
testB <- matrix(NA, n.ind, n.primary)
for (i in 1:nrow(testB)){
  for (j in 1:ncol(testB)){
    if(3 %in% obs[i,j,]) {
      testB[i,j] <- sum(obs[i,j,]==3 | obs[i,j,]==2 | obs[i,j,]==1, na.rm = TRUE)  
    }else{
      testB[i,j] <- 0 
    }
  }
}

seenB <- array(NA, c(n.ind, n.primary, max(n.secondary)))
missedB <- array(NA, c(n.ind, n.primary, max(n.secondary)))

for (i in 1:nrow(testB)){
  for (t in 1:ncol(testB)){
    for (j in 1:n.secondary[t]){
      if(testB[i,t] > 1 & obs[i,t,j] <= 3){seenB[i,t,j] <- 1}
      if(testB[i,t] >= 1 & obs[i,t,j] == 4){missedB[i,t,j] <- 1}
    }
  }
}

yesB <- matrix(NA, n.primary, max(n.secondary))
noB <- matrix(NA, n.primary, max(n.secondary))

for (i in 1:nrow(yesB)){
  for (j in 1:ncol(yesB)){
    yesB[i,j] <- sum(seenB[,i,j], na.rm = TRUE)
    noB[i,j] <- sum(missedB[,i,j], na.rm = TRUE)
  }
}

totalB <- yesB + noB

# Create vector of first encouner occasions
first <- apply(ch,1,findFirst); first[first == "Inf"] <- NA

# Remove individuals that enter on last occasion (should be none)
ch <- subset(ch, first != n.primary)
first <- subset(first, first != n.primary)

  #####-- Get inits and known states --#####

# Use the uncensored CH data to get sensible inits and known states 
# (leveraging information) from observations made after the right-censor date 
# within the anlysis.

n.primary0 <- 31     # number of primary occasions
n.secondary0 <- c(rep(7, 30),3)     # vwector of secondary occasions
n.ind <- nrow(encounter) #number of birds
index0 <- list(1:7, #occasion index
               8:14,
               15:21,
               22:28,
               29:35,
               36:42,
               43:49,
               50:56,
               57:63,
               64:70,
               71:77,
               78:84,
               85:91,
               92:98,
               99:105,
               106:112,
               113:119,
               120:126,
               127:133,
               134:140,
               141:147,
               148:154,
               155:161,
               162:168,
               169:175,
               176:182,
               183:189,
               190:196,
               197:203,
               204:210,
               211:213)

#make 3d array of encounter (inds X primaries X secondaries)
obs0 <- array(NA, dim = c(n.ind,n.primary0, max(n.secondary0))) #init array
for (i in 1:n.primary0){
  obs0[,i,1:n.secondary0[i]] <- encounter0[,index0[[i]]]
}

#make cap history matrix for primary periods only
ch0 <- matrix(NA, n.ind, n.primary0) #init matrix
for (i in 1:n.ind){
  for (t in 1:n.primary0){
    # print(any(obs[i,t,1:n.secondary[t]] != 4))
    ifelse(any(obs0[i,t,1:n.secondary0[t]] != 4), #if observed at all during primary
           #grab the highest state observed (excluding 4s for unobs)
           ch0[i,t] <- max(obs0[i,t,1:n.secondary0[t]][obs0[i,t,1:n.secondary0[t]] != 4]), 
           ch0[i,t] <- 4) #otherwise record a 4 (not seen)
  }
}

#get known z states based on full data
z0 <- knownStateMS(ch0, 4)
z <- z0[, 1:17]
for(i in 1:length(first)){
  z[i,first[i]] <- 1
}

suppressWarnings(z.init <- initStateMS(ch0, first=first, notseen= 4)[,1:17])
z.init[which(!is.na(z))]<-NA #any known Z need to be NA


####---- Finalize Script ----####

# write out .Rdata with all objects neede by JAGS 
save(z.init, first, ch, n.primary, n.secondary, n.ind, yesW, yesM, yesB, 
     totalW, totalM, totalB, z,  evi,
     file = "output/kiwa_dat_JAGS.rdata")

