
library(ggplot2)
library(tidyr)
library(dplyr)

load("./datasets/pred_pow.Rdata"); out1 <- as.data.frame(t(out[[1]]));
load("./datasets/phase3_robust.Rdata")
out2 <- as.data.frame(t(out[[1]]))
out3 <- as.data.frame(t(out[[2]]))
load('./datasets/Supplementary_simulations/phase3_robust_new.Rdata')
out4 <- as.data.frame(t(out[[1]]))

out1 <- gather(out1)
out2 <- gather(out2)
out3 <- gather(out3)
out4 <- gather(out4)

out1$scenario <- 1
out2$scenario <- 2
out3$scenario <- 3
out4$scenario <- 4
data <- rbind(out1,out2,out3,out4)
data$scenario <- as.factor(data$scenario)
data$value <- ifelse(data$value<=0.05,1,0)
I <- which(data$key %in%  c("p.perm","p.adaptive","p.classical","p.fh01","p.lagged","p.rmst"))
data <- data[I,]
data$key <- factor(data$key,levels=c("p.perm","p.lagged","p.fh01","p.adaptive","p.rmst","p.classical"))

est <- aggregate(data$value,by=list(data$scenario, data$key),function(x) mean(x,na.rm=TRUE))
est$Group.2 <- as.factor(est$Group.2)
est$Group.1 <- as.factor(est$Group.1)
est$Group.2 <- relevel(est$Group.2,ref="p.perm")
#levels(est$Group.2) <- c("Permutation","Adaptive","Mantel","Fleming-Harrington","Lagged","RMST")
est$Group.2 <- factor(est$Group.2,levels = rev(levels(est$Group.2)))
levels(est$Group.1) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4")

temp <- c(rep("black",7),"red")
est$color <- as.factor(est$Group.2=="p.perm")
#est$Group.2 <- factor(levels(est$Group.2)[est$Group.2], levels=c("p.classical",
#																						"p.rmst",
#																						"p.adaptive",
#																						"p.fh01",
#																						"p.lagged",
#																						"p.perm"))
p <- ggplot(est) + geom_point(aes(x=x,y=Group.2,col=color),size=5)+facet_wrap(~Group.1,nrow=2)
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

ggsave("./results/figure_sensitivity_new.pdf",width=7,height=7,plot=p,device = "pdf")



