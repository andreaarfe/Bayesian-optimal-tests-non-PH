
library(survival)
library(compiler)
enableJIT(3)
source("./programs/functions_piecewise_exponential.R")
source("./programs/functions_permutation_test.R")
load("./datasets/phase3.Rdata")
set.seed(91343)

# simulation parameters
Nphase2           <- 180
Nphase3           <- 361
tcens             <- 15
marker.prevalence <- 0.5
NREPS             <- 10000
 
# Kaplan-Meier curve in arm 0
KM0    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==0,])

# Kaplan-Meier curve in arm 1
KM1    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==1,])

# Kaplan-Meier curve of censoring times
KMCens <- survfit(Surv(Time,1-Event)~1,data=phase3)

# Simulation parameters - exponential distribution
lambda0 <- 1/exp(survreg(Surv(Time,Event)~1, data=phase3[phase3$Arm==0,], dist="exponential")$coefficients[1])
lambda1 <- 1/exp(survreg(Surv(Time,Event)~1, data=phase3[phase3$Arm==1,], dist="exponential")$coefficients[1])

# Simulation parameters - lognormal distribution
fit0 <- survreg(Surv(Time,Event)~1, data=phase3[phase3$Arm==0,], dist="lognormal")
mu0 <- fit0$coefficients[1]
sigma0 <- fit0$scale
fit1 <- survreg(Surv(Time,Event)~1, data=phase3[phase3$Arm==1,], dist="lognormal")
mu1 <- fit1$coefficients[1]
sigma1 <- fit1$scale

# Function to simulate from a Kaplan-Meier curve
sim.km <- function(n,KM,fumax=15){
  u     <- runif(n)
  ttt   <- quantile(KM,probs = u)$quantile
  ttt   <- ifelse(is.na(ttt) | ttt>fumax, fumax, ttt)
  event <- as.numeric(ttt<=fumax)
  return(list(time=ttt,event=event))
}

# Function to simulate a trial based on KM curves
sim.trial.rob <- function(N,phase3,KM0,KM1,marker.prevalence,KMCens){
  marker  <- rbinom(N,1,marker.prevalence)
  arm   <- sample(phase3$Arm,replace=TRUE,size=N)
  t0    <- sim.km(N,KM0)$time
  t11   <- sim.km(N,KM1)$time
  tcens <- sim.km(N,KMCens)$time
  time    <- ifelse(arm==1 & marker==1,t11,t0)
  event   <- as.numeric(time<=tcens)
  time    <- pmin(time,tcens)
  return(list(marker=marker,arm=arm,time=time,event=event))
}

# Function to simulate a trial based on exponential rates
sim.trial.exp <- function(N,phase3,lambda0,lambda1,marker.prevalence,KMCens){
  marker  <- rbinom(N,1,marker.prevalence)
  arm   <- sample(phase3$Arm,replace=TRUE,size=N)
  t0    <- rexp(N,rate=lambda0)
  t11   <- rexp(N,rate=lambda1)
  tcens <- sim.km(N,KMCens)$time
  time    <- ifelse(arm==1 & marker==1,t11,t0)
  event   <- as.numeric(time<=tcens)
  time    <- pmin(time,tcens)
  return(list(marker=marker,arm=arm,time=time,event=event))
}

# Function to simulate a trial based on log-normal distributions
sim.trial.aft <- function(N,phase3,mu0,sigma0,mu1,sigma1,marker.prevalence,KMCens){
  marker  <- rbinom(N,1,marker.prevalence)
  arm   <- sample(phase3$Arm,replace=TRUE,size=N)
  t0    <- rlnorm(N, meanlog=mu0, sdlog=sigma0)
  t11   <- rlnorm(N, meanlog=mu1, sdlog=sigma1)
  tcens <- sim.km(N,KMCens)$time
  time    <- ifelse(arm==1 & marker==1,t11,t0)
  event   <- as.numeric(time<=tcens)
  time    <- pmin(time,tcens)
  return(list(marker=marker,arm=arm,time=time,event=event))
}

# Function to apply all compared tests
tests <- function(time,event,arm,marker,breaks,
                  # parameters for z=0
                  alpha00,beta00,alpha10,beta10,
                  # parameters for z=1
                  alpha01,beta01,alpha11,beta11){
  # Marker strata 
  I0 <- which(marker==0)
  I1 <- which(marker==1)
  
  # stratified Cox model
  pval.stratcox <- summary(coxph(Surv(time,event)~arm+strata(marker)))$coefficients[5]
  
  # Cox models with bonferroni correction
  pval.cox0 <- summary(coxph(Surv(time[I0],event[I0])~arm[I0]))$coefficients[5]
  pval.cox1 <- summary(coxph(Surv(time[I1],event[I1])~arm[I1]))$coefficients[5]
  pval.bonf.cox <- 2*min(pval.cox0, pval.cox1)
  
  # log-normal AFT models with bonferroni correction
  pval.aft0 <- summary(survreg(Surv(time[I0],event[I0])~arm[I0], dist="lognormal"))$table[2,4]
  pval.aft1 <- summary(survreg(Surv(time[I1],event[I1])~arm[I1], dist="lognormal"))$table[2,4]
  pval.bonf.aft <- 2*min(pval.aft0, pval.aft1)
  
  # Permutation test
  pval.perm <-  perm_test_marker(time,event,arm,marker,breaks,
                                    # parameters for z=0
                                    alpha00,beta00,alpha10,beta10,
                                    # parameters for z=1
                                    alpha01,beta01,alpha11,beta11)$pval
  return(c(pval.stratcox,
           pval.bonf.cox,
           pval.bonf.aft,
           pval.perm))
}

###  Runs the simulation 
pvals           <- matrix(0, nrow=NREPS, ncol=4)
pvals.exp       <- matrix(0, nrow=NREPS, ncol=4)
pvals.aft       <- matrix(0, nrow=NREPS, ncol=4)

colnames(pvals) <-  c("Stratified Cox",
                      "Cox, Bonferroni",
                      "AFT, Bonferroni",
                      "Permutation")
colnames(pvals.exp) <-  c("Stratified Cox",
                      "Cox, Bonferroni",
                      "AFT, Bonferroni",
                      "Permutation")
colnames(pvals.aft) <-  c("Stratified Cox",
                          "Cox, Bonferroni",
                          "AFT, Bonferroni",
                          "Permutation")
for(i in 1:NREPS){
  print(i)
  # Simulates a phase II study
  phase2_dsets <- sim.trial.rob(Nphase2,phase3,KM0,KM1,marker.prevalence,KMCens)
  time   <- phase2_dsets$time
  arm    <- phase2_dsets$arm
  event  <- phase2_dsets$event
  marker <- phase2_dsets$marker

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
  # Scenario 1 (KM)
  phase3_dsets <- sim.trial.rob(Nphase3,phase3,KM0,KM1,marker.prevalence,KMCens)
  time   <- phase3_dsets$time
  arm    <- phase3_dsets$arm
  event  <- phase3_dsets$event
  marker <- phase3_dsets$marker
  pvals[i,]  <- tests(time,event,arm,marker,breaks,
                    # parameters for z=0
                    alpha00,beta00,alpha10,beta10,
                    # parameters for z=1
                    alpha01,beta01,alpha11,beta11)

  # Scenario 2 (exp)
  phase3_dsets <- sim.trial.exp(Nphase3,phase3,lambda0,lambda1,marker.prevalence,KMCens)
  time   <- phase3_dsets$time
  arm    <- phase3_dsets$arm
  event  <- phase3_dsets$event
  marker <- phase3_dsets$marker
  pvals.exp[i,]  <- tests(time,event,arm,marker,breaks,
                          # parameters for z=0
                          alpha00,beta00,alpha10,beta10,
                          # parameters for z=1
                          alpha01,beta01,alpha11,beta11)
  
  # Scenario 3 (AFT)
  phase3_dsets <- sim.trial.aft(Nphase3,phase3,mu0,sigma0,mu1,sigma1,
                                marker.prevalence,KMCens)
  time   <- phase3_dsets$time
  arm    <- phase3_dsets$arm
  event  <- phase3_dsets$event
  marker <- phase3_dsets$marker
  pvals.aft[i,]  <- tests(time,event,arm,marker,breaks,
                          # parameters for z=0
                          alpha00,beta00,alpha10,beta10,
                          # parameters for z=1
                          alpha01,beta01,alpha11,beta11)
}

# Results
sink(file="./results/marker_sims.txt")
print("Scenario 1 - sim from KM")
colMeans(pvals<=0.05)
print("Scenario 2 - sim from exp dist")
colMeans(pvals.exp<=0.05)
print("Scenario 3 - sim from log-normal")
colMeans(pvals.aft<=0.05)
sink() 

save(list=c("pvals", "pvals.aft", "pvals.exp"),
     file="./datasets/phase3_marker_sim.Rdata")




