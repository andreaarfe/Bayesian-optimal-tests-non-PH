
library(survRM2)
library(parallel)
source("./programs/functions_simulations_phase2.R")
source("./programs/functions_weighted_log_rank.R")
source("./programs/functions_permutation_test.R")
source("./programs/functions_piecewise_exponential.R")
source("./programs/functions_simulates_trial_from_KM.R")
load("./datasets/phase3.Rdata")
set.seed(321741)

# sample size of the simulated phase 2 trial
NPHASE2 <- 180

# Number of simulated trials per sample size
NSIM <- 10000

# sample size of simulated phase 3 trials
#NTRIALSAMP <- seq(50,500,50)
NTRIALSAMP <- 361

# Number of random permutation to draw for estimating the permutation test p-value
NPERM <- 1e3

# Censoring time for the simulated trials 
tcens <- 15

# Randomization probability
RAND <- 2/3

# Kaplan-Meier curve in arm 0
KM0    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==0,])

# Kaplan-Meier curve in arm 1
KM1    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==1,])

# Kaplan-Meier curve of censoring times
KMCens <- survfit(Surv(Time,1-Event)~1,data=phase3)


# function to simulate a trial, apply all tests, and return the p-values
sim.results <- function(n,
                        tcens,
                        randprob=2/3,
                        nperm=1e3,
                        phase3data,
                        NPHASE2,
                        KM0,
                        KM1,
                        KMCens){
  
  # generates the phase 2 data
  # New phase 2 study 
  phase2 <- sim.trial.KM(NPHASE2,KM0,KM1,KMCens,phase3data)
  
  # Fix the cutpoints and gets the posterior parameters for the piecewise exponential model
  post <- get.post(times  = phase2$Time,
                   arms   = phase2$Arm,
                   events = phase2$Event)
  breaks       <- post$breaks
  alpha.0.post <- post$alpha.0.post
  alpha.1.post <- post$alpha.1.post
  beta.0.post  <- post$beta.0.post
  beta.1.post  <- post$beta.1.post
  
  # simulates a phase 3 trial
  trial <- sim.trial.phase3(n,tcens,randprob,
                            alpha.0.post,beta.0.post,
                            alpha.1.post,beta.1.post,
                            breaks)
  times <- trial$times
  arm   <- trial$arm
  event <- trial$event
  
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
                           nsim=nperm)$pval
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
               source("./programs/functions_piecewise_exponential.R")})
clusterSetRNGStream(cl1,602314)
out <- vector(mode="list",length=length(NTRIALSAMP))
for(n in seq_along(out)){
  clusterExport(cl1,"n")
  out[[n]] <- parSapply(cl1, 1:NSIM, FUN=function(i) sim.results(NTRIALSAMP[n],tcens,
                                                                 randprob=RAND,
                                                                 nperm=NPERM,
                                                                 phase3,
                                                                 NPHASE2,
                                                                 KM0,
                                                                 KM1,
                                                                 KMCens))
}
stopCluster(cl1)

# save the simulation results and the simulation sample sizes
params <- list(NTRIALSAMP=NTRIALSAMP)
save(list=c("out","params"),file="./datasets/pred_pow.Rdata")



