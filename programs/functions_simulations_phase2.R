
source("./programs/functions_piecewise_exponential.R")

# function to simulate a pair of hazard functions from the posterior given the phase 2 data
sim.haz <- function(alpha.0.post,beta.0.post,
                    alpha.1.post,beta.1.post,
                    breaks){
  lambda <- matrix(0,nrow=2,ncol=length(breaks))
  for(j in seq_along(breaks)){
      lambda[1,j] <- rgamma(1,shape=alpha.0.post[j],rate=beta.0.post[j])
      lambda[2,j] <- rgamma(1,shape=alpha.1.post[j],rate=beta.1.post[j])
  }
  return(lambda)
}

# function to simulate a trial a given pair of hazard function
# censoring is assumed fixed at time tcens in the simulated trials
sim.trial <- function(n,tcens,randprob=0.5,
                      lambda.0,lambda.1,
                      breaks){
  # Simulates the arm indicators
  arm <- rbinom(n,1,randprob)
  # simulates the event times
  evtimes <- vector(mode="numeric",length=n)
  n.0 <- sum(arm==0)
  n.1 <- n - n.0
  evtimes[arm==0] <- sim.piecewise.exp(n.0,lambda.0,breaks)
  evtimes[arm==1] <- sim.piecewise.exp(n.1,lambda.1,breaks)
  # applies censoring
  times <- pmin(evtimes,tcens)
  event <- as.numeric(evtimes<=tcens)
  # outputs the results
  return(list(times=times,event=event,arm=arm))
}

# function to simulate a phase 3 trial from the predictive distribution given the phase 2 data
sim.trial.phase3 <- function(n,tcens,randprob=0.5,
                             alpha.0.post,beta.0.post,
                             alpha.1.post,beta.1.post,
                             breaks){
  # simulates the hazard function from the posterior of the phase 2 data
  haz <- sim.haz(alpha.0.post,beta.0.post,
                 alpha.1.post,beta.1.post,
                 breaks)
  lambda.0 <- haz[1,]
  lambda.1 <- haz[2,]
  # simulates the phase 3 trial
  trial <- sim.trial(n,tcens,randprob,
                     lambda.0,lambda.1,
                     breaks)
  # outputs the simulation results
  return(trial)
}





