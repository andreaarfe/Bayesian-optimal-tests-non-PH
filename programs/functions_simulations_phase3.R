
# function to simulate phase 3 data 
sim.exp.phase3 <- function(NPHASE3  = 361, # same as original phase 3
                       cens     = 15,  # same as in the other simulations
                       randprob = 2/3, # randomization probability
                       tau      = NULL, # length of delayed treatment effects
                       rate.1   = NULL, # rate during induction period
                       rate.2.1 = NULL, # rate after induction period, arm 1
                       rate.2.0 = NULL  # rate after induction period, arm 0
){
  arm <- rbinom(NPHASE3,1,randprob)
  t1 <- rexp(NPHASE3,rate = 1)
  t <- t1
  I <- which(t1<=tau*rate.1)
  t[I] <- t1[I]/rate.1
  I <- which(t1>tau*rate.1 & arm==1)
  t[I] <- tau + (t1[I]-tau*rate.1)/rate.2.1
  I <- which(t1>tau*rate.1 & arm==0)
  t[I] <- tau + (t1[I]-tau*rate.1)/rate.2.0
  tcens <- pmin(t,cens)
  ecens <- as.numeric(t<=cens)
  return(list(Time=tcens,Event=ecens,Arm=arm))
}


