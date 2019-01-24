
library(survival)
set.seed(5445)

### Loads the IPD data used in Alexander et al. NEMJ 2018 (source Lorenzo Trippa)
ipd <- read.csv("./digitized survival curves/IPD_Set_clean.csv",header = TRUE)

# Selecte the digitized IPD for Checkmate 141, OS (Figure 1A of Ferris et al. NEMJ 2016)
# Study arms: nivolumab vs. standard of care
# Randomization was with 2:1 ratio for nivolumab
I <- which(ipd$Trial=="141Checkmate" & ipd$Outcome=="OS" & ipd$Figure=="1A")
data <- ipd[I,c("Time","Event","Arm")]
data$Arm <- ifelse(data$Arm=="SOC",0,1) # 0/1 arm indicator, 1=nivolumab

# Full phase 3 data
phase3 <- data

# Descriptives and Kaplan-Meier curves
pdf("./results/ipd_phase3.pdf")
km3 <- survfit(Surv(Time,Event)~Arm,data=phase3)
capture.output(survdiff(Surv(Time,Event)~Arm,data=phase3)
               ,file="./results/ipd_phase3.txt") 
plot(km3,main="Phase 3")
dev.off()

# Save datasets
save("phase3",file="./datasets/phase3.Rdata")
