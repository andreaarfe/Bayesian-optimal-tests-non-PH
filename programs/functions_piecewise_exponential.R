
# Function to compute the survival function associated of a piecewise exponential distribution
surv.piecewise.exp <- function(times,lambda,breaks){
  delta <- c(diff(breaks),Inf)
  surv <- vector(mode="numeric",length = length(times))
  for(i in seq_along(times)){
    tt <- pmin(delta,pmax(0,times[i]-breaks))
    surv[i] <- exp(-sum(tt*lambda))
  }
  return(surv)
}

# Function to sample from a piecewise exponential distribution
sim.piecewise.exp <- function(n,lambda,breaks){
  m <- length(breaks)
  haz <- lambda[-m]*diff(breaks)
  cumhaz <- cumsum(c(0,haz))
  e <- rexp(n,rate=1)
  I <- findInterval(e,cumhaz)
  times  <- breaks[I]+(e-cumhaz[I])/lambda[I]
  return(times)
}

# Function to compute the breakpoints and posterior parameters for the piecewise exponential model
get.post <- function(times  = NULL,
                     events = NULL,
                     arms   = NULL,
                     breaks = NULL,
                     alpha.0 = rep(0.001,times=length(breaks)),
                     beta.0  = rep(0.001,times=length(breaks)),
                     alpha.1 = rep(0.001,times=length(breaks)),
                     beta.1  = rep(0.001,times=length(breaks))){
  
  # Fix the cutpoints at selected quantiles of the survival distribution, if not provided
  if(is.null(breaks)){
    #KM <- survival::survfit(Surv(times,events)~1)
    #breaks <- quantile(KM,probs=seq(0,1,1/10))$quantile
    #breaks <- breaks[!is.na(breaks)]  
    breaks <-  c(0,quantile(times[arms==0],probs=(1:4)/5))
  }
  
  # Events per interval
  I.0 <- findInterval(times[events==1 & arms==0],
                      breaks,all.inside = FALSE,left.open = TRUE, rightmost.closed = TRUE)
  I.1 <- findInterval(times[events==1 & arms==1],
                      breaks,all.inside = FALSE,left.open = TRUE, rightmost.closed = TRUE)
  I.0 <- factor(I.0,levels = 1:length(breaks))
  I.1 <- factor(I.1,levels = 1:length(breaks))
  ev.0 <- table(I.0)
  ev.1 <- table(I.1)
  
  # Person-time per interval
  delta <- c(diff(breaks),Inf)
  pt.0 <- rep(0,times=length(breaks))
  pt.1 <- rep(0,times=length(breaks))
  if(sum(arms==0)>0) pt.0 <- rowSums(sapply(times[arms==0],function(t) pmin(delta,pmax(0,t-breaks))))
  if(sum(arms==1)>0) pt.1 <- rowSums(sapply(times[arms==1],function(t) pmin(delta,pmax(0,t-breaks))))
  pt.0[is.na(pt.0)] <- 0
  pt.1[is.na(pt.1)] <- 0
  
  # Posterior parameters
  alpha.0.post <- alpha.0 + ev.0
  beta.0.post  <- beta.0  + pt.0
  alpha.1.post <- alpha.1 + ev.1
  beta.1.post  <- beta.1  + pt.1
  
  # Returns the output
  return(list(breaks=breaks,
              alpha.0.post=alpha.0.post,
              alpha.1.post=alpha.1.post,
              beta.0.post=beta.0.post,
              beta.1.post=beta.1.post))
}
  

