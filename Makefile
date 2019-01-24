
DATPATH=./datasets
FUNPATH=./programs
RESPATH=./results

all:           data figures simulations

data:          $(DATPATH)/phase3.Rdata
               
simulations:   $(DATPATH)/pred_pow.Rdata \
               $(DATPATH)/phase3_boot_power.Rdata \
               $(DATPATH)/phase3_robust.Rdata \
               $(RESPATH)/marker_sims.txt
               
figures:       $(RESPATH)/figure_phase2_pred_power.pdf \
               $(RESPATH)/figure_sensitivity.pdf \
               $(RESPATH)/figure_phase3_post_plot.pdf

clean:
	rm -fv $(DATPATH)/*.Rdata
	rm -fv $(RESPATH)/*.*
	rm -fv $(FUNPATH)/*.Rout
	rm -fv ./*.Rout
	rm -fv ./nohup.out

cleanlogs:
	rm -fv $(FUNPATH)/*.Rout
	rm -fv ./*.Rout
	rm -fv ./nohup.out


####### Datasets


$(DATPATH)/phase3.Rdata: $(FUNPATH)/data_prepIPD.R
	R CMD BATCH $(FUNPATH)/data_prepIPD.R


####### Simulations


# Predictive power
$(DATPATH)/pred_pow.Rdata: $(FUNPATH)/simulation_pred_pow.R \
			   $(FUNPATH)/functions_permutation_test.R \
			   $(FUNPATH)/functions_simulations_phase2.R \
			   $(FUNPATH)/functions_piecewise_exponential.R \
			   $(FUNPATH)/functions_weighted_log_rank.R \
			   $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/simulation_pred_pow.R

# Bootstrap power estimates
$(DATPATH)/phase3_boot_power.Rdata: \
                           $(FUNPATH)/simulation_bootstrap_power.R \
			   $(FUNPATH)/functions_permutation_test.R \
			   $(FUNPATH)/functions_simulations_phase2.R \
			   $(FUNPATH)/functions_piecewise_exponential.R \
			   $(FUNPATH)/functions_weighted_log_rank.R \
			   $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/simulation_bootstrap_power.R

# Robustness analysis
$(DATPATH)/phase3_robust.Rdata: \
                           $(FUNPATH)/simulation_power_robustness_exp.R \
			   $(FUNPATH)/functions_permutation_test.R \
		           $(FUNPATH)/functions_simulations_phase2.R \
			   $(FUNPATH)/functions_piecewise_exponential.R \
			   $(FUNPATH)/functions_weighted_log_rank.R \
			   $(FUNPATH)/functions_simulations_phase3.R \
			   $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/simulation_power_robustness_exp.R

$(RESPATH)/marker_sims.txt: $(FUNPATH)/simulation_marker_KM.R \
			    $(FUNPATH)/functions_permutation_test.R \
                            $(FUNPATH)/functions_piecewise_exponential.R \
                            $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/simulation_marker_KM.R

####### Figures

$(RESPATH)/figure_phase3_post_plot.pdf: $(FUNPATH)/figure_posterior_plot.R \
					$(FUNPATH)/functions_piecewise_exponential.R 
					$(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/figure_posterior_plot.R

$(RESPATH)/figure_phase2_pred_power.pdf: $(FUNPATH)/figure_phase2_pred_power.R \
                                         $(DATPATH)/phase3_boot_power.Rdata 
	R CMD BATCH $(FUNPATH)/figure_phase2_pred_power.R

	
$(RESPATH)/figure_sensitivity.pdf: $(DATPATH)/pred_pow.Rdata \
                                   $(DATPATH)/phase3_robust.Rdata \
				   $(FUNPATH)/figure_sensitivity_analysis.R
	R CMD BATCH $(FUNPATH)/figure_sensitivity_analysis.R
     
     
