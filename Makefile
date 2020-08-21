
DATPATH=./datasets
FUNPATH=./programs
RESPATH=./results

all:           data figures simulations sim_phase2_size sim_hypeparams

data:          $(DATPATH)/phase3.Rdata
               
simulations:   $(DATPATH)/pred_pow.Rdata \
               $(DATPATH)/phase3_boot_power.Rdata \
               $(DATPATH)/phase3_robust.Rdata \
               $(RESPATH)/marker_sims.txt \
               $(DATPATH)/Supplementary_simulations/phase3_marker_sim_2.Rdata \
               $(RESPATH)/Supplementary_results/marker_sims_2.txt \
               $(DATPATH)/Supplementary_simulations/phase3_marker_sim_2.Rdata \
               $(DATPATH)/Supplementary_simulations/phase3_robust_new.Rdata
               
figures:       $(RESPATH)/figure_phase2_pred_power.pdf \
               $(RESPATH)/figure_sensitivity_new.pdf \
               $(RESPATH)/figure_phase3_post_plot.pdf \
               $(RESPATH)/figure_marker.pdf \
               $(RESPATH)/Supplementary_results/figure_marker_2.pdf \
               $(RESPATH)/Supplementary_results/figure_surv_funct_sims.pdf \
               $(RESPATH)/Supplementary_results/figure_phaseII_size.pdf
							 
sim_phase2_size:  $(DATPATH)/Supplementary_simulations/sim_phaseII_size_40.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_60.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_80.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_100.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_120.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_140.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_160.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_180.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_200.Rdata \
									$(DATPATH)/Supplementary_simulations/sim_phaseII_size_220.Rdata 

sim_hypeparams: $(DATPATH)/Supplementary_simulations/sim_hyper_u1_v1.Rdata \
                $(DATPATH)/Supplementary_simulations/sim_hyper_u10_v10.Rdata


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

DEPENDENCIES=$(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R \
											     $(FUNPATH)/functions_permutation_test.R \
												   $(FUNPATH)/functions_simulations_phase2.R \
												   $(FUNPATH)/functions_piecewise_exponential.R \
												   $(FUNPATH)/functions_weighted_log_rank.R \
												   $(DATPATH)/phase3.Rdata
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
												   $(FUNPATH)/functions_simulates_trial_from_KM.R \
												   $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/simulation_power_robustness_exp.R

$(DATPATH)/phase3_marker_sim.Rdata \
$(RESPATH)/marker_sims.txt: $(FUNPATH)/simulation_marker_KM.R \
														$(FUNPATH)/functions_permutation_test.R \
                            $(FUNPATH)/functions_piecewise_exponential.R \
                            $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/simulation_marker_KM.R

####### Figures

$(RESPATH)/figure_phase3_post_plot.pdf: $(FUNPATH)/figure_posterior_plot.R \
																				$(FUNPATH)/functions_piecewise_exponential.R \
																				$(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/figure_posterior_plot.R

$(RESPATH)/figure_phase2_pred_power.pdf: $(FUNPATH)/figure_phase2_pred_power.R \
																				 $(DATPATH)/phase3_boot_power.Rdata 
	R CMD BATCH $(FUNPATH)/figure_phase2_pred_power.R

	
# $(RESPATH)/figure_sensitivity.pdf: $(DATPATH)/pred_pow.Rdata \
#																	 $(DATPATH)/phase3_robust.Rdata \
#																	 $(FUNPATH)/figure_sensitivity_analysis.R
#	R CMD BATCH $(FUNPATH)/figure_sensitivity_analysis.R

$(RESPATH)/figure_sensitivity_new.pdf: \
																	 $(DATPATH)/Supplementary_simulations/phase3_robust_new.Rdata \
																	 $(DATPATH)/pred_pow.Rdata \
																	 $(DATPATH)/phase3_robust.Rdata \
																	 $(FUNPATH)/figure_sensitivity_analysis_new.R
	R CMD BATCH $(FUNPATH)/figure_sensitivity_analysis_new.R
     
$(RESPATH)/figure_marker.pdf: $(DATPATH)/phase3_marker_sim.Rdata \
															$(FUNPATH)/figure_marker_sim.R
	R CMD BATCH $(FUNPATH)/figure_marker_sim.R


####### Supplementary analyses and results

# Plot of Kaplan-Meier curse corresponding to the scenarios in Section 7.2
$(RESPATH)/Supplementary_results/figure_surv_funct_sims.pdf: \
													 $(FUNPATH)/Supplementary_analyses/Figure_survival_curves_sec7.2.R \
											     $(FUNPATH)/functions_permutation_test.R \
												   $(FUNPATH)/functions_simulations_phase2.R \
												   $(FUNPATH)/functions_piecewise_exponential.R \
												   $(FUNPATH)/functions_weighted_log_rank.R \
												   $(FUNPATH)/functions_simulations_phase3.R \
												   $(FUNPATH)/functions_simulates_trial_from_KM.R \
												   $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/Supplementary_analyses/Figure_survival_curves_sec7.2.R

# Simulations of marker-stratified design in high-power setting
$(DATPATH)/Supplementary_simulations/phase3_marker_sim_2.Rdata \
$(RESPATH)/Supplementary_results/marker_sims_2.txt: \
														$(FUNPATH)/Supplementary_analyses/simulation_marker_KM_2.R \
														$(FUNPATH)/functions_permutation_test.R \
                            $(FUNPATH)/functions_piecewise_exponential.R \
                            $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/Supplementary_analyses/simulation_marker_KM_2.R


$(RESPATH)/Supplementary_results/figure_marker_2.pdf: \
															$(DATPATH)/Supplementary_simulations/phase3_marker_sim_2.Rdata \
															$(FUNPATH)/Supplementary_analyses/figure_marker_sim_2.R
	R CMD BATCH $(FUNPATH)/Supplementary_analyses/figure_marker_sim_2.R
	
# Simulations varying the phase II sample size
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_40.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=40' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_60.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=60' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_80.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=80' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_100.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=100' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_120.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=120' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_140.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=140' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_160.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=160' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_180.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=180' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_200.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=200' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(DATPATH)/Supplementary_simulations/sim_phaseII_size_220.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args NPHASE2=220' $(FUNPATH)/Supplementary_analyses/simulation_phaseII_samp_size.R 
	
$(RESPATH)/Supplementary_results/figure_phaseII_size.pdf: \
								$(DATPATH)/Supplementary_simulations/sim_phaseII_size_*.Rdata \
								$(FUNPATH)/Supplementary_analyses/figure_phaseII_size.R 
	R CMD BATCH $(FUNPATH)/Supplementary_analyses/figure_phaseII_size.R 

# Simulations varying the gamma hyper-params
$(DATPATH)/Supplementary_simulations/sim_hyper_u1_v1.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args u=1 v=1' $(FUNPATH)/Supplementary_analyses/simulation_hyper_params.R 

$(DATPATH)/Supplementary_simulations/sim_hyper_u10_v10.Rdata: $(DEPENDENCIES)
	R CMD BATCH '--args u=10 v=10' $(FUNPATH)/Supplementary_analyses/simulation_hyper_params.R 
	
# Additional scenarios for robustness analysis of Section 7.2

$(DATPATH)/Supplementary_simulations/phase3_robust_new.Rdata: \
                           $(FUNPATH)/Supplementary_analyses/simulation_power_robustness_exp_2.R \
											     $(FUNPATH)/functions_permutation_test.R \
												   $(FUNPATH)/functions_simulations_phase2.R \
												   $(FUNPATH)/functions_piecewise_exponential.R \
												   $(FUNPATH)/functions_weighted_log_rank.R \
												   $(FUNPATH)/functions_simulations_phase3.R \
												   $(FUNPATH)/functions_simulates_trial_from_KM.R \
												   $(DATPATH)/phase3.Rdata
	R CMD BATCH $(FUNPATH)/Supplementary_analyses/simulation_power_robustness_exp_2.R
