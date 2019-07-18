
library(ggplot2)
library(tidyverse)
load("./datasets/phase3_marker_sim.Rdata")

pvals <- as.data.frame(pvals)
pvals.exp <- as.data.frame(pvals.exp)
pvals.aft <- as.data.frame(pvals.aft)
pvals$scen <- "Scenario 1"
pvals.exp$scen <- "Scenario 2"
pvals.aft$scen <- "Scenario 3"
d <- gather(rbind(pvals, pvals.exp, pvals.aft),
						"Stratified Cox", "Cox, Bonferroni", 
						"AFT, Bonferroni", "Permutation",
						key="test", value="pval")
d$scen <- as.factor(d$scen)
d$test <- as.factor(d$test)
d$color <- d$test=="Permutation"
levels(d$test) <- c("Stratum-specific AFT", "Stratum-specific log-rank",
										"Permutation", "Stratified Cox")

p <- d %>% 
	group_by(scen,test,color) %>%
	summarise(pow=mean(pval<=0.05)) %>%
	arrange(pow, .by_group=TRUE) %>%
	ungroup() %>%
	mutate(test=fct_reorder(test, pow)) %>%
	ggplot(aes(y=test,x=pow,col=color)) +
	facet_wrap(~scen, nrow=1) + 
	geom_point(size=5) + 
	theme_bw() + 
	ylab("Test type") + 
	xlab("Probability of rejection") + 
	scale_color_manual(values=c("black","red")) + 
	theme(text = element_text(size=25),
				legend.position="none",
				panel.spacing=unit(1.5,"line"),
				plot.margin=unit(c(5.5, 10.5, 5.5, 5.5), "points")) +
	scale_x_continuous(breaks = c(0.00, 0.20,0.40,0.6), limits=c(0,0.6))

ggsave("./results/figure_marker.pdf",width=10,height=4,plot=p,device = "pdf")
