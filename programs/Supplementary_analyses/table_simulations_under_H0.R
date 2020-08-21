
library(xtable)
library(tidyverse)

load('./datasets/Supplementary_simulations/sim_no_trt_eff.Rdata')
sims <- as.data.frame(t(out)) %>% gather(key=pval)
tab <- sims %>%
	group_by(pval) %>%
	summarise(alpha = mean(value<=0.05),
						N = n()) %>%
	mutate(SE = sqrt(alpha*(1-alpha)/N),
				 LCL95 = alpha - 1.96*SE,
				 UCL95 = alpha + 1.96*SE) %>%
	select(pval, alpha, LCL95, UCL95) %>%
	xtable(digits=3)
print(tab, include.rownames=FALSE)

