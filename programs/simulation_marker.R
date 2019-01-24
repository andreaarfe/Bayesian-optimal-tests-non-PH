
library(survival)
source("./programs/functions_piecewise_exponential.R")
source("./programs/functions_permutation_test.R")
load("./datasets/phase3.Rdata")
set.seed(79135)

# simulation parameters
Nphase2           <- 500
Nphase3           <- 1000
tcens             <- 17
rand              <- 2/3
marker.prevalence <- 0.5
NREPS             <- 10000

# Event rates
pt1 <- sum(phase3$Time[phase3$Arm  == 1])
ev1 <- sum(phase3$Event[phase3$Arm == 1])
pt0 <- sum(phase3$Time[phase3$Arm  == 0])
ev0 <- sum(phase3$Event[phase3$Arm == 0])
r1 <- ev1/pt1
r0 <- ev0/pt0

#  Runs the simulation 
pval.perm     <- vector(mode="numeric",length = NREPS)
pval.stratcox <- vector(mode="numeric",length = NREPS)
for(i in 1:NREPS){
  print(i)
  # Simulates a phase II study
  marker  <- rbinom(Nphase2,1,marker.prevalence)
  arm     <- rbinom(Nphase2,1,rand)
  time    <- ifelse(arm==1,
                    rexp(Nphase2,rate=marker*r1+(1-marker)*r0),
                    rexp(Nphase2,rate=r0))
  event   <- as.numeric(time<=tcens)
  time    <- pmin(time,tcens)
  
  # Defines the cut-points of the model
  breaks <- c(0,quantile(time[arm==0],probs=(1:4)/5))
  
  # Fit the piecewise exponential model with binary marker
  I0 <- which(marker==0)
  I1 <- which(marker==1)
  data0 <- get.post(time[I0],event[I0],arm[I0],breaks)
  data1 <- get.post(time[I1],event[I1],arm[I1],breaks)
  alpha00 <- data0$alpha.0.post
  beta00  <- data0$beta.0.post
  alpha10 <- data0$alpha.1.post
  beta10  <- data0$beta.1.post
  alpha01 <- data1$alpha.0.post
  beta01  <- data1$beta.0.post
  alpha11 <- data1$alpha.1.post
  beta11  <- data1$beta.1.post
  
  # Simulates a phase III study
  marker  <- rbinom(Nphase3,1,marker.prevalence)
  arm     <- rbinom(Nphase3,1,rand)
  time    <- ifelse(arm==1,
                    rexp(Nphase3,rate=marker*r1+(1-marker)*r0),
                    rexp(Nphase3,rate=r0))
  event   <- as.numeric(time<=tcens)
  time    <- pmin(time,tcens)
  
  # stratified Cox model
  pval.stratcox[i] <- summary(coxph(Surv(time,event)~arm+strata(marker)))$coefficients[5]
  
  # Permutation test
  pval.perm[i] <-  perm_test_marker(time,event,arm,marker,breaks,
                                    # parameters for z=0
                                    alpha00,beta00,alpha10,beta10,
                                    # parameters for z=1
                                    alpha01,beta01,alpha11,beta11)$pval
}

# Results
pow <- c(mean(pval.perm<=0.05),mean(pval.stratcox<=0.05))
names(pow) <- c("Permutation","Stratified Cox")
sink(file="./results/marker_sims.txt")
print(pow)
sink()
