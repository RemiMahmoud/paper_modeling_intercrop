#  Modeling intercrop systems

This repo contains the code related to the paper submitted on ...

Its is structured in multiple folders.

## R

Contains the different R scripts related to different parts of the analysis.

In this README, " $\rightarrow$ " means that the script produces the elements (datasets / tables) after the $\rightarrow$

### 1-field_experiments_description

*0-explore_modeling_possibilities* : Identify experiments and traits suited for the analysis

*0-describe_fields_experiments* : Quick desc of the experiments involved in the Models $\rightarrow$ doc/table_environment.tex


### 2-modeling

*0-Get_control_managements_IC_SC* : Find control management in sole crop for each intercrop experimental unit (same N, same cultivar etc.) $\rightarrow$ data/management_control_IC_SC.rds

*1-Compute_LER_and_biodiversity_effect* :  Compute Land Equivalent Ratio and biodiversity decomposition for each experimental unit $\rightarrow$ data/data_BE_CE_SE_LER.rds

*1-Compute_NNI* : Compute the Nitrogen Nutrition Index for each experimental unit $\rightarrow$ data/data_NNI.rds

*1-Non_parametric_fit_height_and_biomass* : Use smoothing splines to reduce the dimension of growth curves (height / biomass) $\rightarrow$ data/data_(biomass|height|cover)_fit.rds  

*2-Find_maximal_value_for_each_variable* : Get the maximum biomass and height value $\rightarrow$ data/data_asymp.rds  

*3-Get_values_LAI_SLA_Cover* : Get the maximum LAI / SLA and cover values $\rightarrow$ data/data_SLA_LAI_cover.rds

*4-Impute_values_and_compute_differences* : Carry out the imputation of missing data and compute differences $\rightarrow$ data/data_all_diff_with_imputations.rds

*5-Run_all_models* :  Run the Mixed effect random forest models for yields and (P)LER $\rightarrow$ models stored in models/MERF/(wt_fababean|wt_pea|all_species)/*.rds

*5-Run_all_models_CV* : Run the Mixed effect random forest models for yields and (P)LER, with cross validation $\rightarrow$ models stored in models/MERF/models_train_test_cv/(wt_fababean|wt_pea|all_species)/*.rds

*5-Run_all_models_N_seed* :  Run the Mixed effect random forest models for N content in seeds (cereal) $\rightarrow$ models stored in models/N_seed/(wt_fababean|wt_pea|all_species)/*.rds

*5-Run_all_models_N_seed_CV* : Run the Mixed effect random forest models for N content in seeds (cereal)  , with cross validation $\rightarrow$ models stored in models/N_seed_CV/(wt_fababean|wt_pea|all_species)/*.rds

*6-Figure_correlation_between_performance_variables* :  Figure IV.1 Manuscript : correlation between performance variables $\rightarrow$ figures/plot_correlations_performance.png

### 3-analysis_models

*6_Analysis_models_cv_LongituRF* : Predicted vs obs plots for yield/(P)LER models with cross validation $\rightarrow$ figures/MERF/models_cv*

*6_Analysis_models_LongituRF* : Fitted vs obs plots for yield/(P)LER models ; and importance plots for each variable $\rightarrow$ figures/MERF/plot_fitted_vs_obs_all_models.png

*6_Analysis_models_N_seed_LongituRF* : Fitted vs obs plots for N seed models ; and importance plots for each variable $\rightarrow$ figures/N_seed/plot_fitted_vs_obs_all_models.png

*7-Marginal_plot_LongituRF* : Plot and save all marginal plots for the yield/(P)LER models $\rightarrow$ figures/MERF/marginal_plots*

*7-Marginal_plot_N_seed_LongituRF* : Plot and save all marginal plots for the Nseed models $\rightarrow$ figures/N_seed/marginal_plots*

*8-Summary_models* : Plot the summary of the models as a tile plot, with variable ranked by their importance, and signed by the $\tau$ with the response variable $\rightarrow$ figures/MERF/
figures_summaries_models*

*9-Figure_error_models* : Plot the error (RMSE / EF or any chosen metric) $\rightarrow$ figures/MERF/error_models_cv_all.png

*10-Figure_random_effects_model* : Caterpillar plot of the random effects $\rightarrow$ figures/MERF/plot_RE_models.png ; doc/table_variance_imputation.tex

*11-Code_for_further_analysis* : Miscellaneous to look closer at relationships

*12-Figures_PhD_defence* : Big figures for defence phd

*13-Figures_Distribution_number_selections* : Figure of how many times a variables was selected xx times for each model


## Data

Contains differents datasets used during the different steps of the analysis

## Models

Contains the different fitted models, in .rds files (can be read with `readr::read_rds(*models/model_xxx.rds*)`)

## Doc

*exploration_covariates* : Distribution, description of the covariates

*boruta_NNI_illustration* : Figure of binomial distribution to explain Boruta and figure of critical N curve
