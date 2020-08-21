
library(ggplot2)
library(tidyr)
library(dplyr)
source("./programs/functions_simulations_phase2.R")
source("./programs/functions_weighted_log_rank.R")
source("./programs/functions_permutation_test.R")
source("./programs/functions_piecewise_exponential.R")
source("./programs/functions_simulations_phase3.R")
source("./programs/functions_simulates_trial_from_KM.R")
load('./datasets/Supplementary_simulations/phase3_robust_new.Rdata')
out4 <- as.data.frame(t(out[[2]]))
out4 <- gather(out4)

times <- seq(0,15,0.5)
rate.2.1 <- simparams[[2]]$rate.2.1
rate.2.0 <- simparams[[2]]$rate.2.0
rate.1.0 <- simparams[[2]]$rate.1.0
rate.1.1 <- simparams[[2]]$rate.1.1
S0 <- surv.piecewise.exp(times, lambda = c(rate.1.0, rate.2.0), breaks = c(0, 2))
S1 <- surv.piecewise.exp(times, lambda = c(rate.2.0, rate.2.1), breaks = c(0, 2))
pdf("./results/Supplementary_results/crossing_KM.pdf",
		height=6, width=6)
plot(times, S0, type='l', ylim=c(0,1), 
		 col='blue', lwd=3,
		 xlab = 'Months of follow-up',
		 ylab = 'Survival probability',
		 main='Supplementary scenario', 
		 cex.lab=1.5, cex.main=2, cex.axis=2)
points(times, S1, type='l', ylim=c(0,1), 
			 col='orange', lwd=3)
points(times, S0, type='l', ylim=c(0,1), 
			 col='blue', lwd=3)
dev.off()

data <- out4
data$value <- ifelse(data$value<=0.05,1,0)
I <- which(data$key %in%  c("p.perm","p.adaptive","p.classical","p.fh01","p.lagged","p.rmst"))
data <- data[I,]
data$key <- factor(data$key,levels=c("p.perm","p.lagged","p.fh01","p.adaptive","p.rmst","p.classical"))

est <- aggregate(data$value,by=list(data$key),function(x) mean(x,na.rm=TRUE))
est$Group.1 <- as.factor(est$Group.1)
est$Group.1 <- relevel(est$Group.1,ref="p.perm")
est$Group.1 <- factor(est$Group.1,levels = rev(levels(est$Group.1)))

temp <- c(rep("black",7),"red")
est$color <- as.factor(est$Group.1=="p.perm")
#est$Group.2 <- factor(levels(est$Group.2)[est$Group.2], levels=c("p.classical",
#																						"p.rmst",
#																						"p.adaptive",
#																						"p.fh01",
#																						"p.lagged",
#																						"p.perm"))
p <- ggplot(est) + geom_point(aes(x=x,y=Group.1,col=color),size=5)
p <- p + scale_color_manual(values=c("black","red"))
p <- p + theme_bw()
p <- p + xlab("Probability of rejection") + ylab("Test type")
p <- p + theme(text = element_text(size=25),
							 panel.spacing=unit(1.5,"line"),
               legend.position="none",
               axis.text.y = element_text(colour = temp))
#p <- p + xlim(0,1)
p <- p + scale_x_continuous(breaks = c(0.20,0.80),minor_breaks = c(0,0.4,0.6,1), limits=c(0,1))
p <- p + scale_y_discrete(labels=rev(expression("Permutation",
																								"Lagged (10%)",
																								G^"0,1",
                                                "Adaptive",
                                                "RMST",
																								"Mantel")) )

ggsave("./results/Supplementary_results/figure_sims_crossing.pdf",width=6,height=6,plot=p,device = "pdf")



