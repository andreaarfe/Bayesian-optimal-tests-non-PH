
library(survRM2)
library(parallel)
source("./programs/functions_simulations_phase2.R")
source("./programs/functions_weighted_log_rank.R")
source("./programs/functions_permutation_test.R")
source("./programs/functions_piecewise_exponential.R")
source("./programs/functions_simulations_phase3.R")
source("./programs/functions_simulates_trial_from_KM.R")
load("./datasets/phase3.Rdata")
set.seed(3419047)

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

# Grid of time points for the plots
times <- seq(0,15,0.5)

### Phase 3 simulation parameters for scenarios 2 and 3

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

# Scenario 3: same rate until 8 months of follow-up, then different rates
ev.3.1   <- sum(phase3$Event[phase3$Time<=8])
pt.3.1   <- sum(pmin(phase3$Time,8))
ev.3.2.1 <- sum(phase3$Event[phase3$Arm==1 & phase3$Time>8])
pt.3.2.1 <- sum(pmax(phase3$Time[phase3$Arm==1]-8,0))
ev.3.2.0 <- sum(phase3$Event[phase3$Arm==0 & phase3$Time>8])
pt.3.2.0 <- sum(pmax(phase3$Time[phase3$Arm==0]-8,0))
rate.3.1 <- ev.3.1/pt.3.1
rate.3.2.1 <- ev.3.2.1/pt.3.2.1
rate.3.2.0 <- ev.3.2.0/pt.3.2.0

# Scenario 4: no delays (proportional hazards)
ev.4.1   <- sum(phase3$Event)
pt.4.1   <- sum(phase3$Time)
ev.4.2.1 <- sum(phase3$Event[phase3$Arm==1])
pt.4.2.1 <- sum(phase3$Time[phase3$Arm==1])
ev.4.2.0 <- sum(phase3$Event[phase3$Arm==0])
pt.4.2.0 <- sum(phase3$Time[phase3$Arm==0])
rate.4.1 <- ev.4.1/pt.4.1
rate.4.2.1 <- ev.4.2.1/pt.4.2.1
rate.4.2.0 <- ev.4.2.0/pt.4.2.0

### Sample survival curves for scenario 1

# function to simulate a trial, apply all tests, and return the p-values
sim.results <- function(n,
												tcens,
												randprob=2/3,
												nperm=1e3,
												phase3data,
												NPHASE2,
												KM0,
												KM1,
												KMCens,
												times){
	
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
	
	# Survival probabilities
	lambda0 <- matrix(0,nrow=1,ncol=length(breaks))
	lambda1 <- matrix(0,nrow=1,ncol=length(breaks))
	for(j in seq_along(breaks)){
		lambda0[1,j] <- rgamma(1,shape=alpha.0.post[j],rate=beta.0.post[j])
		lambda1[1,j] <- rgamma(1,shape=alpha.1.post[j],rate=beta.1.post[j])
	}
	out <- vector(mode='list', length = 2)
	out[[1]] <- surv.piecewise.exp(times = times, lambda = lambda0, breaks = breaks)
	out[[2]] <- surv.piecewise.exp(times = times, lambda = lambda1, breaks = breaks)
	return(out)
}
S <- lapply(1:NSIM, function(i) sim.results(NTRIALSAMP[n],tcens,
																								randprob=RAND,
																								nperm=NPERM,
																								phase3,
																								NPHASE2,
																								KM0,
																								KM1,
																								KMCens,
																			 times))

#uclS0 <- apply(S0samp, 1, quantile, probs=0.975)
#lclS0 <- apply(S0samp, 1, quantile, probs=0.025)
#uclS1 <- apply(S1samp, 1, quantile, probs=0.975)
#lclS1 <- apply(S1samp, 1, quantile, probs=0.025)

pdf(file = './results/Supplementary_results/figure_surv_funct_sims.pdf',
		width = 9, height = 9)
par(mfrow=c(2,2))
S0samp <- sapply(seq_along(S), function(i) S[[i]][[1]])
S1samp <- sapply(seq_along(S), function(i) S[[i]][[2]])
S0 <- apply(S0samp, 1, mean)
S1 <- apply(S1samp, 1, mean)
plot(times, S0, type='l', ylim=c(0,1), 
		 col='blue', lwd=3,
		 xlab = 'Months of follow-up',
		 ylab = 'Survival probability',
		 main='Scenario 1 (predictive distributions)')
#polygon(c(times, rev(times)), c(uclS0,rev(lclS0)), 
#				col=rgb(0,0,0.5, 0.2), border=NA)
#polygon(c(times, rev(times)), c(uclS1,rev(lclS1)), 
#				col=rgb(0.5,0.4,0, 0.2), border=NA)
points(times, S1, type='l', ylim=c(0,1), 
			 col='orange', lwd=3)
points(times, S0, type='l', ylim=c(0,1), 
			 col='blue', lwd=3)

S0 <- surv.piecewise.exp(times, lambda = c(rate.2.1, rate.2.2.0), breaks = c(0, 2))
S1 <- surv.piecewise.exp(times, lambda = c(rate.2.1, rate.2.2.1), breaks = c(0, 2))
plot(times, S0, type='l', ylim=c(0,1), 
		 col='blue', lwd=3,
		 xlab = 'Months of follow-up',
		 ylab = 'Survival probability',
		 main='Scenario 2 (2 months delay)')
points(times, S1, type='l', ylim=c(0,1), 
		 col='orange', lwd=3)
points(times, S0, type='l', ylim=c(0,1), 
			 col='blue', lwd=3)

S0 <- surv.piecewise.exp(times, lambda = c(rate.3.1, rate.3.2.0), breaks = c(0, 8))
S1 <- surv.piecewise.exp(times, lambda = c(rate.3.1, rate.3.2.1), breaks = c(0, 8))
plot(times, S0, type='l', ylim=c(0,1), 
		 col='blue', lwd=3,
		 xlab = 'Months of follow-up',
		 ylab = 'Survival probability',
		 main = 'Scenario 3 (8 months delay)')
points(times, S1, type='l', ylim=c(0,1), 
			 col='orange', lwd=3)
points(times, S0, type='l', ylim=c(0,1), 
			 col='blue', lwd=3)

S0 <- surv.piecewise.exp(times, lambda = c(rate.4.1, rate.4.2.0), breaks = c(0, 0))
S1 <- surv.piecewise.exp(times, lambda = c(rate.4.1, rate.4.2.1), breaks = c(0, 0))
plot(times, S0, type='l', ylim=c(0,1), 
		 col='blue', lwd=3,
		 xlab = 'Months of follow-up',
		 ylab = 'Survival probability',
		 main = 'Scenario 4 (no delays, proportional hazards)')
points(times, S1, type='l', ylim=c(0,1), 
			 col='orange', lwd=3)
points(times, S0, type='l', ylim=c(0,1), 
			 col='blue', lwd=3)
dev.off()
