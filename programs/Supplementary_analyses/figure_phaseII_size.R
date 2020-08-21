
library(tidyverse)

sizes <- seq(40,220,20)
envs <- vector(mode='list', length = length(sizes))
tests <- vector(mode='list', length = length(sizes))
for(i in seq_along(sizes)){
	envs[[i]] <- new.env()
	fname <- paste0('./datasets/Supplementary_simulations/sim_phaseII_size_',
									sizes[i],
									'.Rdata')
	load(fname, envs[[i]])
	res <- mean(ifelse(envs[[i]]$out$sims <= 0.05, 1, 0))
	tests[[i]] <- data.frame(NPHASE2=envs[[i]]$out$NPHASE2, 
													 pow=res)
}
allres <- bind_rows(tests)

p <- ggplot(allres, aes(x=NPHASE2, y=pow)) +
	#geom_ribbon(aes(ymin=allres$pow - sqrt(allres$pow*(1-allres$pow)/10000)*1.96, 
	#								ymax=allres$pow + sqrt(allres$pow*(1-allres$pow)/10000)*1.96), 
	#						alpha=0.3) +
	geom_line() +
	geom_point() +
	theme_bw() +
	ggtitle('Power estimates for different phase II sample size') +
	xlab('Phase II sample size') +
	ylab('Power') +
	ylim(0.8,0.95) +
	theme(text = element_text(size=10)) +
	scale_x_continuous(breaks=sizes)

ggsave(plot=p, 
			 filename = './results/Supplementary_results/figure_phaseII_size.pdf',
			 width = 5,
			 height = 4)
