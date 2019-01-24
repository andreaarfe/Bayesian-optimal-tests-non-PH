
library(survival)
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

# Function to simulate from a Kaplan-Meier curve
sim.km <- function(n,KM,fumax=15){
  u     <- runif(n)
  ttt   <- quantile(KM,probs = u)$quantile
  ttt   <- ifelse(is.na(ttt) | ttt>fumax, fumax, ttt)
  event <- as.numeric(ttt==fumax)
  return(list(time=ttt,event=event))
}

# Function to simulate a trial
sim.trial.rob <- function(N,phase3,KM0,KM1,KMCens,marker.prevalence,rand){
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

#  Runs the simulation 
pval.perm     <- vector(mode="numeric",length = NREPS)
pval.perm.no  <- vector(mode="numeric",length = NREPS)
pval.stratcox <- vector(mode="numeric",length = NREPS)
pval.cox0     <- vector(mode="numeric",length = NREPS)
pval.cox1     <- vector(mode="numeric",length = NREPS)
pval.bonf     <- vector(mode="numeric",length = NREPS)
phase2_dsets  <- vector(mode="list",length = NREPS)
phase3_dsets  <- vector(mode="list",length = NREPS)
for(i in 1:NREPS){
  #print(i)
  # Simulates a phase II study
  phase2_dsets[[i]] <- sim.trial.rob(Nphase2,phase3,KM0,KM1,KMCens,marker.prevalence,rand)
  time   <- phase2_dsets[[i]]$time
  arm    <- phase2_dsets[[i]]$arm
  event  <- phase2_dsets[[i]]$event
  marker <- phase2_dsets[[i]]$marker
  
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
  phase3_dsets[[i]] <- sim.trial.rob(Nphase3,phase3,KM0,KM1,KMCens,marker.prevalence,rand)
  time   <- phase3_dsets[[i]]$time
  arm    <- phase3_dsets[[i]]$arm
  event  <- phase3_dsets[[i]]$event
  marker <- phase3_dsets[[i]]$marker
  
  # stratified Cox model
  pval.stratcox[i] <- summary(coxph(Surv(time,event)~arm+strata(marker)))$coefficients[5]
  
  # Cox-models with bonferroni correction
  I0 <- which(marker==0)
  I1 <- which(marker==1)
  pval.cox0[i] <- summary(coxph(Surv(time[I0],event[I0])~arm[I0]))$coefficients[5]
  pval.cox1[i] <- summary(coxph(Surv(time[I1],event[I1])~arm[I1]))$coefficients[5]
  
  # Permutation test
  pval.perm[i] <-  perm_test_marker(time,event,arm,marker,breaks,
                                    # parameters for z=0
                                    alpha00,beta00,alpha10,beta10,
                                    # parameters for z=1
                                    alpha01,beta01,alpha11,beta11)$pval
}
pval.bonf <- 2*pmin(pval.cox0,pval.cox1)
#pval.sidak <- 1-(1-pmin(pval.cox0,pval.cox1))^2

# Results
pow <- c(mean(pval.perm<=0.05), 
        	mean(pval.stratcox<=0.05), 
	          mean(pval.bonf<=0.05)
                #,mean(pval.sidak<=0.05)
         )
names(pow) <- c("Permutation",
                "Stratified Cox",
                "Bonferroni"
                #,"Sidak"
                )

save(list="pow",file="./datasets/phase3_marker_sim.Rdata")




