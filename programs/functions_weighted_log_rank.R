
library(survival)
library(YPmodel)

# Function to perform the weighted two-samples log-rank test
# length(weight) == number of the distinct ordered event times across both groups
weighted.logrank <- function(time, event, group, weights=1){
  event1 <- event*group
  event0 <- event*(1-group)
  tt <- sort(unique(time))
  exits0 <- sapply(tt, function(t) sum((group==0)[time==t]))
  exits1 <- sapply(tt, function(t) sum((group==1)[time==t]))
  n1 <- c(sum(group==1), sum(group==1)-cumsum(exits1)); n1 <- n1[-length(n1)]
  n0 <- c(sum(group==0), sum(group==0)-cumsum(exits0)); n0 <- n0[-length(n0)]
  d0 <- sapply(tt,function(t) sum(event0[time==t]))
  d1 <- sapply(tt,function(t) sum(event1[time==t]))
  evtimes <- which(d0+d1>0)
  n0_ev <- n0[evtimes]
  n1_ev <- n1[evtimes]
  d0_ev <- d0[evtimes]
  d1_ev <- d1[evtimes]
  d <- d0_ev+d1_ev
  n <- n0_ev+n1_ev
  obs <- d1_ev
  exp <- d*n1_ev/n
  num <- exp-obs
  denum <- ((n-d)*d*n0_ev*n1_ev)/((n-1)*n^2); denum[n==1] <- 0
  X2 <- sum(weights*num)^2/sum((weights^2)*denum)
  p <- pchisq(X2,1,lower.tail = FALSE)
  out <- list("chi-square"=X2,"p-value"=p)
  return(out)
}

# Lagged log-rank test
lagged.logrank <- function(time,event,group,lag=0){
  tt <- sort(unique(time[event==1]))
  weights <- as.numeric(tt>=lag)
  weighted.logrank(time,event,group,weights)
}

# Fleming-Harrington log-rank test
flemharr.logrank <- function(time,event,group,rho=0,delta=0){
  S <- survfit(Surv(time,event)~1)
  S_ev <- S$surv[S$n.event>0]
  weights <- (S_ev^rho)*((1-S_ev)^delta)
  weighted.logrank(time,event,group,weights)
}

# adaptive log-rank of Yang and Prentice (Biometrics 2010)
adaptive.logrank <- function(time,event,group){
  test <- suppressWarnings(YPmodel.adlgrk(data.frame(V1=time,V2=event,V3=group)))
  return(test$pval)
}

# Piecewise-weighted log-rank test
# piecewise.logrank <- function(time,event,group,cuts,weights){
#   tt <- unique(time[event==1])
#   I <- findInterval(tt,cuts)
#   w <- weights[I+1]
#   weighted.logrank(time,event,group,w)
# }

