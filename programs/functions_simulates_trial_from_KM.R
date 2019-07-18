
library(survival)

sim.trial.KM <- function(n,           # Sample size of the simulated trial
                         KM0,         # Kaplan-Meier curve in arm 0
                         KM1,         # Kaplan-Meier curve in arm 1
                         KMCensoring, # Kaplan-Meier curve of censoring times
                         phase3data,   # Original phase 3 dataset
                         fumax=15    # Maximum duration of the simulated trial
                                      # 15 = max follow-up in SOC arm of CheckMate 141
                         ){
  
  # Generates patients
  I <- sample.int(length(phase3data$Arm),size=n,replace=TRUE)
  Arm <- phase3data$Arm[I]
  n1 <- sum(Arm) 
  n0 <- n - n1
  Time  <- rep(0,times=n)
  Event <- rep(0,times=n)
  
  # Generates censoring times
  ucens    <- runif(n)
  cens <- quantile(KMCens,probs = ucens)$quantile
  cens <- ifelse(is.na(cens) | cens>fumax, fumax, cens)
  
  # Generates follow-up times in arm 0
  if(n0>0){
    u0 <- runif(n0)
    t0 <- quantile(KM0,probs = u0)$quantile
    Time[Arm==0]  <- ifelse(is.na(t0) | t0 > cens[Arm==0], cens[Arm==0], t0)
    Event[Arm==0] <- ifelse(is.na(t0) | t0 > cens[Arm==0], 0, 1)
  }
  
  # Generates follow-up times in arm 1
  if(n1>0){
    u1 <- runif(n1)
    t1 <- quantile(KM1,probs = u1)$quantile
    Time[Arm==1]  <- ifelse(is.na(t1) | t1 > cens[Arm==1], cens[Arm==1], t1)
    Event[Arm==1] <- ifelse(is.na(t1) | t1 > cens[Arm==1], 0, 1)
  }
  
  # Simulated trial 
  return(list(Arm=Arm,Event=Event,Time=Time))
}

# #Tests the function
# load("./datasets/phase3.Rdata")
# #Kaplan-Meier curve in arm 0
# KM0    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==0,])
# #Kaplan-Meier curve in arm 1
# KM1    <- survfit(Surv(Time,Event)~1,data=phase3[phase3$Arm==1,])
# #Kaplan-Meier curve of censoring times
# KMCens <- survfit(Surv(Time,1-Event)~1,data=phase3)
# #Simulates a trial and compares the KM curves with the original data
# trial <- sim.trial.KM(1e5,KM0,KM1,KMCens,phase3)
# plot(survfit(Surv(Time,Event)~Arm,data=trial),xlim=c(0,15))
# plot(survfit(Surv(Time,Event)~Arm,data=phase3),xlim=c(0,15))


