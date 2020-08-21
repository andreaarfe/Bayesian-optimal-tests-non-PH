
library(tidyverse)

load('./datasets/Supplementary_simulations/sim_hyper_u1_v1.Rdata')
d1 <- out$sims
load('./datasets/Supplementary_simulations/sim_hyper_u10_v10.Rdata')
d10 <- out$sims
load('./datasets/phase3_boot_power.Rdata')
d10em3 <- out[6,]
allres <- data.frame(pow=c(mean(as.numeric(d1<=0.05)),
													 mean(as.numeric(d10<=0.05)),
													 mean(as.numeric(d10em3<=0.05))),
										 pvar=c(1, 0.1, 1000))
allres$pvar <- as.factor(allres$pvar)

p <- ggplot(allres, aes(x=pvar, y=pow)) +
	#geom_ribbon(aes(ymin=allres$pow - sqrt(allres$pow*(1-allres$pow)/10000)*1.96, 
	#								ymax=allres$pow + sqrt(allres$pow*(1-allres$pow)/10000)*1.96), 
	#						alpha=0.3) +
	geom_point() +
	theme_bw() +
	ggtitle('Power estimates for the tailored permutation test') +
	xlab('Prior variance') +
	ylab('Power') +
	ylim(0.85,0.95) +
	theme(text = element_text(size=10)) 

ggsave(plot=p, 
			 filename = './results/Supplementary_results/figure_hyperparams.pdf',
			 width = 5,
			 height = 4)
