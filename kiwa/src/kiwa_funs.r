####################################
####                            ####
####       KIWA Functions       ####
####      Scott Yanco, PhD      ####
####    scott.yanco@yale.edu    ####
####                            ####
####################################


####---- Description ----####

# This script contains custom functions used in preparing data for analysis.
# Currently all fucntions in this sccript are sourced only by the `prep_data.r`
# script.

####---- notin ----####

# Simple negation of `in` function
"%notin%" <- Negate("%in%")


####---- findFirst ----####

# Returns which elements in a vector != 4 (the not seen code) and then returns 
# the minimum of those values
findFirst <- function(x) {
  min(which(x != 4))
}


####---- knownStateMS ----####

#create known observed states for z from data
knownStateMS <- function(ch, notseen){
  state <- ch
  for (i in 1:dim(ch)[1]){
    #winter
    if(1 %in% ch[i,]) {
      W1 <- min(which(ch[i,]==1))
      W2 <- max(which(ch[i,]==1))
      state[i,W1:W2] <- 1
      state[i,W1] <- NA
    }
    
    #migration
    if(2 %in% ch[i,]){
      M1 <- min(which(ch[i,]==2))
      M2 <- max(which(ch[i,]==2))
      state[i,M1:M2] <- 2
    }
    
    #breeding
    if(3 %in% ch[i,]){
      B1 <- min(which(ch[i,]==3))
      B2 <- max(which(ch[i,]==3))
      state[i,B1:B2] <- 3
      
    }
  }
  state[state==notseen] <- NA
  return(state)
}


####---- initStateMS ----####
initStateMS <- function(ch, notseen, first){
  #init empty known state matrix
  ks <- matrix(NA, nrow = dim(ch)[1], ncol = dim(ch)[2])
  ch[ch==notseen]<-NA
  for(i in 1:dim(ch)[1]){ #loop thorugh ind
    for(j in (first[i]+1):dim(ch)[2]){ #and time steps, starting at first cap
      if(is.na(ch[i,j])){ #if ind not observed that occasion...
        #assign it the highest state observed thus far
        before <- max(na.omit(ch[i,1:j]))
        after <- min(na.omit(ch[i,j:dim(ch)[2]]))
        ks[i,j] <- floor(mean(c(before,after)))
        if(is.infinite(ks[i,j])){ks[i,j]<-7} #if not observed after, guess dead
      } #if
    } #j
  } #i
  return(ks)
} # fxn
