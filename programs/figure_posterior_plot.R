
library(survival)
library(ggplot2)
source("./programs/functions_piecewise_exponential.R")
load("./datasets/phase3.Rdata")
set.seed(3241)

# Fix the cutpoints at selected quantiles of the survival distribution
KM <- survfit(Surv(Time,Event)~1,data=phase3)

# Posterior parameters
post <- get.post(times=phase3$Time,
                 events=phase3$Event,
                 arms=phase3$Arm)
alpha.0.post <- post$alpha.0.post
beta.0.post  <- post$beta.0.post
alpha.1.post <- post$alpha.1.post
beta.1.post  <- post$beta.1.post
breaks       <- post$breaks

# samples the posterior distribution
nsamp <- 100
tt <- seq(0,max(phase3$Time),0.1)
lambda.0 <- matrix(0,nrow=nsamp,ncol=length(breaks))
lambda.1 <- matrix(0,nrow=nsamp,ncol=length(breaks))
S.0 <- matrix(0,nrow=nsamp,ncol=length(tt))
S.1 <- matrix(0,nrow=nsamp,ncol=length(tt))
for(i in 1:nsamp){
  for(j in seq_along(breaks)){
    lambda.0[i,j] <- rgamma(1,shape=alpha.0.post[j],rate=beta.0.post[j])
    lambda.1[i,j] <- rgamma(1,shape=alpha.1.post[j],rate=beta.1.post[j])
  }
  S.0[i,] <- surv.piecewise.exp(tt,lambda.0[i,],breaks)
  S.1[i,] <- surv.piecewise.exp(tt,lambda.1[i,],breaks)
}

# plot of the posterior samples
# pdf("./results/figure_phase3_post_plot.pdf")
# par(mar=c(4.1,4.4,1,1))
# plot(survfit(Surv(Time,Event)~Arm,data=phase3),
#      lwd=3,
#      xlab="Months of follow-up",
#      ylab="Survival time",
#      cex.lab = 2,
#      cex.axis=1.5)
# points(tt,colMeans(S.0),type="l",col=rgb(0,0,1),lwd=3)
# points(tt,colMeans(S.1),type="l",col=rgb(1,0.4,0),lwd=3)
# for(i in 1:nsamp){
#   points(tt,S.0[i,],type="l",col=rgb(0,0,1,alpha=0.05),lwd=2)
#   points(tt,S.1[i,],type="l",col=rgb(1,0.4,0,alpha=0.05),lwd=2)
# }
#text(7,0.1,labels=c("Standard of care"),col=rgb(0,0,1),cex=1.5)
#text(14,0.5,labels=c("Nivolumab"),col=rgb(1,0.4,0),cex=1.5)
#dev.off()


est <- data.frame(S=c(colMeans(S.0),colMeans(S.1)),
                  time=c(tt,tt),
                  group=as.factor(rep(c(0,1),each=length(tt))))
fit <- survfit(Surv(Time,Event)~Arm,data=phase3)
d <- data.frame(time=fit$time,surv=fit$surv)
d$group <- as.factor(c(rep(1,times=fit$strata[1]),rep(0,times=fit$strata[2])))
d <- rbind(d,data.frame(time=c(0,0),surv=c(1,1),group=c(0,1)))
p <- ggplot() + geom_step(data=d,aes(x=time,y=surv,group=group),size=1.1) + ylim(0,1)
p <- p + geom_line(data=est,aes(x=time,y=S,col=group),size=1.1)
p <- p + theme_classic()
p <- p + xlab("Months of follow-up") + ylab("Survival probability")
p <- p + theme(text = element_text(size=25))
p <- p + scale_color_manual(values=c(rgb(0,0,1),rgb(1,0.4,0))) 
p <- p + guides(color=FALSE)
p <- p + annotate("text",x=14,y=0.5,label="Nivolumab",col=rgb(1,0.4,0),size=7)
p <- p + annotate("text",x=7,y=0.1,label="Standard of care",col=rgb(0,0,1),size=7)
p <- p + geom_vline(xintercept = breaks[-1], linetype = "dotted")
ggsave("./results/figure_phase3_post_plot.pdf",height = 4,width = 7,device = "pdf", plot=p)



