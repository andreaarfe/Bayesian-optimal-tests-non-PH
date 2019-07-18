
library(ggplot2)
load("./datasets/phase3_boot_power.Rdata"); 
#NSAMP <- params$NTRIALSAMP

# prepares the power estimates for plotting
#pow <- data.frame()
#for(i in seq_along(NSAMP)){
  #pvals <- out[[i]]
  #powest <- rowMeans(pvals<=0.05)
  #newrows <- cbind(N=rep(NSAMP[i],times=length(powest)),
  #                 pow=powest,
  #                 type=names(powest))
  #pow <- rbind(pow,newrows)
#}
pow <- data.frame(pow=colMeans(t(out<=0.05)))
pow$type <- rownames(pow)
rownames(pow) <- NULL
I <- which(pow$type %in% c("p.perm",
                           "p.adaptive",
                           "p.classical",
                           "p.fh01",
                           "p.lagged",
                           "p.rmst"))
pow <- pow[I,]
pow$type <- factor(pow$type,levels=c("p.perm",
                                     "p.adaptive",
                                     "p.classical",
                                     "p.fh01",
                                     "p.lagged",
                                     "p.rmst"))
#pow$N <- as.numeric(levels(pow$N))[pow$N]
#pow$pow <- as.numeric(levels(pow$pow))[pow$pow]
pow$type <- relevel(pow$type,ref="p.perm")
#levels(pow$type) <- c("Permutation","Adaptive","Mantel","FH(0,1)","FH(1,0)","FH(1,1)","Lagged","RMST")
levels(pow$type) <-  c("Permutation",
                       "Adaptive",
                       "Mantel",
                       "FH (rho=1)",
                       "Lagged (10%)",
                       "RMST")
pow$type <- factor(pow$type, levels=rev(levels(pow$type)))


# plots the power estimates
#pdf("./results/figure_phase2_pred_power.pdf",width=7,height=7)
# p <- ggplot(pow,aes(N,pow,col=type,shape=type))+ylim(0.25,1)+scale_shape_manual(values=1:nlevels(pow$type)) 
# p <- p + geom_line(size=1.1)+geom_point(size=5)
# p <- p + xlab("Phase III Sample size") + ylab("Bayesian Expected Power")
# p <- p + theme_classic()
# p <- p + theme(text = element_text(size=25))
# p <- p + theme(legend.position = c(0.81, 0.25))
# p <- p + guides(shape=guide_legend(title="Test type:",keywidth=2,keyheight = 1.5),
#                 col=guide_legend(title="Test type:")) 
#plot(p)

pow$color <- as.factor(pow$type=="Permutation")
temp <- c(rep("black",7),"red")
pow <- pow[rev(order(pow$pow)),]
pow$type <- factor(pow$type, levels=rev(levels(pow$type)[pow$type]))

p <- ggplot(pow) + geom_point(aes(x=pow,y=type,color=color),size=4)
p <- p + scale_color_manual(values=c("black","red"))
p <- p + theme_bw()
p <- p + xlab("Probability of rejection") + ylab("Test type")
p <- p + theme(text = element_text(size=25),
               legend.position="none",
               axis.text.y = element_text(colour = temp))
#p <- p + scale_x_continuous(breaks = c(0.00,0.20,0.40,0.60,0.80,1.00),limits=c(0,1))
p <- p + scale_x_continuous(breaks = c(0.00,0.20,0.40,0.50,0.60,0.70,0.80,0.90,1.00),limits=c(0.5,1))
p <- p + scale_y_discrete(labels=rev(expression("Permutation",
                                                "Lagged (10%)",
                                                G^"0,1",
                                            "Adaptive",
                                            "RMST",
                                            "Mantel"
                                            )) )

ggsave("./results/figure_phase2_pred_power.pdf",width=7,height=4,plot=p,device="pdf")
#dev.off()



