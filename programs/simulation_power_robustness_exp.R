
library(survRM2)
library(parallel)
source("./programs/functions_simulations_phase2.R")
source("./programs/functions_weighted_log_rank.R")
source("./programs/functions_permutation_test.R")
source("./programs/functions_piecewise_exponential.R")
source("./programs/functions_simulations_phase3.R")
source("./programs/functions_simulates_trial_from_KM.R")
load("./datasets/phase3.Rdata")
set.seed(325235)

# Phase 2 sample size
NPHASE2 <- 180 

# Number of replicates per scenario
NSIM <- 10000

# Kaplan-Meier curve in arm 0
KM0    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==0,])

# Kaplan-Meier curve in arm 1
KM1    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==1,])

# Kaplan-Meier curve of censoring times
KMCens <- survfit(Surv(Time,1-Event)~1,data=phase3)

# Phase 3 simulation parameters
#simparams <- vector(mode="list",length=5)
simparams <- vector(mode="list",length=2)

# Scenario 1: different rates (no delay)
# ev.1.1 <- sum(phase3$Event[phase3$Arm==1])
# ev.1.0 <- sum(phase3$Event[phase3$Arm==0])
# pt.1.1 <- sum(phase3$Time[phase3$Arm==1])
# pt.1.0 <- sum(phase3$Time[phase3$Arm==0])
# rate.1.1 <- ev.1.1/pt.1.1
# rate.1.0 <- ev.1.0/pt.1.0
# simparams[[1]] <- list(tau = 0,
#                        rate.1 = 1,
#                        rate.2.1 = rate.1.1,
#                        rate.2.0 = rate.1.0)

# Scenario 2: same rate until 2 months of follow-up, then different rates
ev.2.1   <- sum(phase3$Event[phase3$Time<=2])
pt.2.1   <- sum(pmin(phase3$Time,2))
ev.2.2.1 <- sum(phase3$Event[phase3$Arm==1 & phase3$Time>2])
pt.2.2.1 <- sum(pmax(phase3$Time[phase3$Arm==1]-2,0))
ev.2.2.0 <- sum(phase3$Event[phase3$Arm==0 & phase3$Time>2])
pt.2.2.0 <- sum(pmax(phase3$Time[phase3$Arm==0]-2,0))
rate.2.1 <- ev.2.1/pt.2.1
rate.2.2.1 <- ev.2.2.1/pt.2.2.1
rate.2.2.0 <- ev.2.2.0/pt.2.2.0
simparams[[1]] <- list(tau = 2,
                       rate.1 = rate.2.1,
                       rate.2.1 = rate.2.2.1,
                       rate.2.0 = rate.2.2.0)

# Scenario 3: same rate until 4 months of follow-up, then different rates
# ev.3.1   <- sum(phase3$Event[phase3$Time<=4])
# pt.3.1   <- sum(pmin(phase3$Time,4))
# ev.3.2.1 <- sum(phase3$Event[phase3$Arm==1 & phase3$Time>4])
# pt.3.2.1 <- sum(pmax(phase3$Time[phase3$Arm==1]-4,0))
# ev.3.2.0 <- sum(phase3$Event[phase3$Arm==0 & phase3$Time>4])
# pt.3.2.0 <- sum(pmax(phase3$Time[phase3$Arm==0]-4,0))
# rate.3.1 <- ev.3.1/pt.3.1
# rate.3.2.1 <- ev.3.2.1/pt.3.2.1
# rate.3.2.0 <- ev.3.2.0/pt.3.2.0
# simparams[[3]] <- list(tau = 4,
#                        rate.1 = rate.3.1,
#                        rate.2.1 = rate.3.2.1,
#                        rate.2.0 = rate.3.2.0)

# Scenario 4: same rate until 8 months of follow-up, then different rates
ev.4.1   <- sum(phase3$Event[phase3$Time<=8])
pt.4.1   <- sum(pmin(phase3$Time,8))
ev.4.2.1 <- sum(phase3$Event[phase3$Arm==1 & phase3$Time>8])
pt.4.2.1 <- sum(pmax(phase3$Time[phase3$Arm==1]-8,0))
ev.4.2.0 <- sum(phase3$Event[phase3$Arm==0 & phase3$Time>8])
pt.4.2.0 <- sum(pmax(phase3$Time[phase3$Arm==0]-8,0))
rate.4.1 <- ev.4.1/pt.4.1
rate.4.2.1 <- ev.4.2.1/pt.4.2.1
rate.4.2.0 <- ev.4.2.0/pt.4.2.0
simparams[[2]] <- list(tau = 8,
                       rate.1 = rate.4.1,
                       rate.2.1 = rate.4.2.1,
                       rate.2.0 = rate.4.2.0)

# Scenario 5: same rates (no effect)
# ev.5 <- sum(phase3$Event)
# pt.5 <- sum(phase3$Time)
# rate.5 <- ev.5/pt.5
# simparams[[5]] <- list(tau = 100,
#                        rate.1 = rate.5,
#                        rate.2.1 = 1,
#                        rate.2.0 = 1)

# Plots of KM curves simulated from each scenario
# NPHASE3 <- 100000
# for(i in 1:5){
#   outfile <- paste0("./results/km_scenario",i,".pdf")
#   pdf(file = outfile)
#   ph3 <- sim.exp.phase3(NPHASE3 = NPHASE3,
#                     tau = simparams[[i]]$tau,
#                     rate.1 = simparams[[i]]$rate.1,
#                     rate.2.1 = simparams[[i]]$rate.2.1,
#                     rate.2.0 = simparams[[i]]$rate.2.0)
#   plot(survfit(Surv(Time,Event)~Arm,data=ph3),main=paste0("Scenario ",i),
#        xlab="Follow-up duration (months)",
#        ylab="Survival probability")
#   KM <- survfit(Surv(Time,Event)~Arm,data=phase3)
#   points(c(0,KM$time[1:KM$strata[1]]),c(1,KM$surv[1:KM$strata[1]]),
#          type="s",col="red",lty=2)
#   points(c(0,KM$time[(1+KM$strata[1]):(length(KM$time))]),c(1,KM$surv[(1+KM$strata[1]):length(KM$time)]),
#          type="s",col="blue",lty=2)
#   dev.off()
# }

# function to simulate a phase 3 trial, apply all tests, and return the p-values
sim.results <- function(scenario,
                        simparams,
                        NPERM=1e3,
                        phase3data,
                        NPHASE2,
                        KM0,
                        KM1,
                        KMCens){
  # New phase 2 study 
  phase2 <- sim.trial.KM(NPHASE2,KM0,KM1,KMCens,phase3data)

  # breaks and posterior parameters
  post <- get.post(times  = phase2$Time,
                   arms   = phase2$Arm,
                   events = phase2$Event)
  breaks       <- post$breaks
  alpha.0.post <- post$alpha.0.post
  alpha.1.post <- post$alpha.1.post
  beta.0.post  <- post$beta.0.post
  beta.1.post  <- post$beta.1.post

  # simulated phase 3 trial 
  trial <- sim.exp.phase3(NPHASE3  = 361,
                          tau      = simparams[[scenario]]$tau,
                          rate.1   = simparams[[scenario]]$rate.1,
                          rate.2.1 = simparams[[scenario]]$rate.2.1,
                          rate.2.0 = simparams[[scenario]]$rate.2.0)
  times <- trial$Time
  arm   <- trial$Arm
  event <- trial$Event

  # Lag-time for the lagged-log-rank (excludes the first 10% of follow-up times)
  lag <- quantile(times[arm==0],probs=0.10)
  
  # cut-point for RMST (lowest maximum follow-up time between arms) 
  tau <-  min(max(times[arm==0]),max(times[arm==1]))
  
  # applies the tests
  p.classical <- weighted.logrank(times,event,arm)$"p-value"
  p.lagged    <- lagged.logrank(times,event,arm,lag=lag)$"p-value"
  p.fh01      <- flemharr.logrank(times,event,arm,rho=0,delta=1)$"p-value"
  #p.fh10      <- flemharr.logrank(times,event,arm,rho=1,delta=0)$"p-value"
  #p.fh11      <- flemharr.logrank(times,event,arm,rho=1,delta=1)$"p-value"
  p.adaptive  <- adaptive.logrank(times,event,arm)
  p.rmst      <- rmst2(times,event,arm,tau=tau)$unadjusted.result[1,"p"]
  #p.rmst      <- rmst2(times,event,arm)$unadjusted.result[1,"p"]
  p.perm      <- perm_test(times,event,arm,breaks,
                           alpha.0.post,beta.0.post,
                           alpha.1.post,beta.1.post,
                           nsim=NPERM)$pval
  # outputs the results
  out <- c(p.classical,
           p.lagged,
           p.fh01,
           #p.fh10,
           #p.fh11,
           p.adaptive,
           p.rmst,
           p.perm)
  names(out) <- c("p.classical",
                  "p.lagged",
                  "p.fh01",
                  #"p.fh10",
                  #"p.fh11",
                  "p.adaptive",
                  "p.rmst",
                  "p.perm")
  return(out)
}


# run the simulation
NCLUST <- 4 # number of clusters on which to run the simulation
cl1<-makeCluster(NCLUST)
clusterExport(cl1, ls())
clusterEvalQ(cl1,
             {library(survRM2)
               library(parallel)
               source("./programs/functions_simulations_phase2.R")
               source("./programs/functions_weighted_log_rank.R")
               source("./programs/functions_permutation_test.R")
               source("./programs/functions_piecewise_exponential.R")
               source("./programs/functions_simulations_phase3.R")
               source("./programs/functions_simulates_trial_from_KM.R")
             })
clusterSetRNGStream(cl1,3419047)
out <- vector(mode="list",length=length(simparams))
for(n in seq_along(out)){
  clusterExport(cl1,"n")
  out[[n]] <- parSapply(cl1, 1:NSIM, FUN=function(i) sim.results(n,
                                                                 simparams,
                                                                 NPERM=1e3,
                                                                 phase3,
                                                                 NPHASE2,
                                                                 KM0,
                                                                 KM1,
                                                                 KMCens))
}
stopCluster(cl1)

# save the simulation results and the simulation sample sizes
save(list=c("out","simparams"),file="./datasets/phase3_robust.Rdata")


