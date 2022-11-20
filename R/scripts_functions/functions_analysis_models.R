
# Functions to get fitted values of the GMERF models
fun_fitted_gmerf <- function(gm, data, y="yield_cereal"){
  
  fitted_values <- gm %>% fitted.gmerf(type = "mu")
  data_gmerf_fit_vs_obs <- tibble(observed = data[[y]], fitted = fitted_values)
  
  return(data_gmerf_fit_vs_obs)
}


fun_get_fitted_observed <-function(data_models, y ){
  data_models %>% 
    mutate(fitted = map2(mod.gmerf, data, fun_fitted_gmerf, y= y )) %>% 
    unnest(fitted) %>% 
    select(.imp, fitted, observed) %>%
    arrange(.imp, observed) %>% 
    group_by(.imp) %>% 
    mutate(id_obs = 1:n()) %>% 
    ungroup %>% 
    return
}


# Function to get the importance if the random forests
fun_compute_importance_gmerf <- function(gm){
  
  # Get a tibble for the importance of the variables
  imp <- gm$forest.model$importance
  
  tibble_importance_gmerf  <- tibble(variable = rownames(imp), importance = as.vector(imp)) 
  return(tibble_importance_gmerf)
  
}

# Function to summarise the importance of the different covariates in the models 
get_importance <- function(data_models){
  
  
  sign_importance_variable_numeric <-  data_models %>% 
    unnest(data) %>% 
    select(c(contains("diff")), contains("NNI"), contains("cover"), response = any_of(matches("yield|LER|BE|CE|SE", ignore.case = FALSE))) %>%
    pivot_longer(-response, names_to = "variable") %>% 
    group_by(variable) %>% 
    summarise(cor_kendal = cor(response,value, method = "kendall"), .groups = "drop") %>%
    mutate(sign_cor = fct_recode(factor(sign(cor_kendal)), "NEG" = "-1", "POS"  = "1"))
  
  
  data_models %>% 
    mutate(importance = map(mod.gmerf, fun_compute_importance_gmerf)) %>% 
    unnest(importance) %>% 
    group_by(variable) %>% 
    summarise(mean_importance = mean(importance), sd_importance = sd(importance), min_importance = mean(importance), max_importance = max(importance),  n = n(), .groups = "drop") %>% 
    left_join(sign_importance_variable_numeric) %>% 
    mutate(variable = stringr::str_replace(variable, "cereal" , "C")) %>%  
    mutate(variable = stringr::str_replace(variable, "legume" , "L")) %>%  
    mutate(variable = forcats::fct_reorder(variable, mean_importance)) %>% 
    return
}



# Fucntion to decompose the variance of the residuals into unexplained and variance due to imputations
variance_decomposition <- function(data_fitted_observed){
  
  summary_fitted_observed <- data_fitted_observed %>%
    group_by(id_obs) %>% 
    summarise(mean_fitted = mean(fitted), observed = unique(observed), .groups = "drop")
  
  SSRT <- sum((data_fitted_observed$fitted - data_fitted_observed$observed)^2)
  
  data_variance_imputation <- data_fitted_observed %>% left_join(summary_fitted_observed) 
  
  M <- length(unique(data_fitted_observed$.imp))
  
  variance_imputation <- sum((data_variance_imputation$mean_fitted - data_variance_imputation$fitted)^2)/SSRT
  variance_unexplained <- M*sum((summary_fitted_observed$mean_fitted - summary_fitted_observed$observed)^2)/SSRT
  
  return(list(SSRT = SSRT, prop_var_imputation = variance_imputation, prop_var_unexplained  = variance_unexplained))
}

# Function to copute the error of the GMERF models
get_error <- function(data_models, list_metrics = c("bias", "RMSE", "RRMSE", "EF"), y = "yield_cereal"){
  
  data_models %>% 
    mutate(fitted = map2(mod.gmerf, data, fun_fitted_gmerf, y=y )) %>% 
    mutate(error.gmerf = map(fitted, evaluate_error)) %>% 
    unnest(error.gmerf) %>% 
    summarise(across(any_of(list_metrics), 
                     .fns = list(mean = mean, sd = sd),
                     .names = "{.col}_{.fn}"), .groups = "drop")
  
}


#####################################
#### Plotting functions #############
#####################################

plot_fitted_vs_observed <- function(data_fitted_observed, summary_fitted_observed, my_title = NULL, lab_x = "observed", lab_y = "fitted", size_text = 18, add_error_label = TRUE, list_error_metrics = c("RMSE", "RRMSE", "EF")){
  
  plot_output <- data_fitted_observed %>% 
    ggplot(aes(x = observed, y= fitted)) + 
    geom_point(alpha=.2) + 
    geom_point(data= summary_fitted_observed, aes(y = mean_fitted), shape = 4) + 
    geom_y_x() + 
    labs(title = my_title, x = lab_x, y = lab_y)+
    theme(title = element_text(size = size_text),
          axis.text = element_text(size = size_text +2))
  
  if(add_error_label){
    
    tibble_error <- data_fitted_observed %>%
      group_by(.imp) %>%
      nest %>%
      ungroup %>%  
      mutate(error.gmerf = map(data, evaluate_error)) %>% 
      unnest(error.gmerf) %>% 
      summarise(across(any_of(list_error_metrics), 
                       .fns = list(mean = mean),
                       .names = "{.col}_{.fn}"), .groups = "drop") %>% 
      round(digits = 2)
    
    
    # label_error <- TeX(paste(c("$\\bar{RMSE}$ = ", tibble_error$RMSE_mean, "\newline", 
    #                        "$\\bar{RRMSE}$ = " , 100*tibble_error$RRMSE_mean, "% \newline",
    #                        "$\\bar{EF}$ = ", 100*tibble_error$EF_mean, "% \newline", 
    #                        "%Var imputation = ", round(variance_decomposition(data_fitted_observed)$prop_var_imputation, digits =2 )*100, "%"), 
    #                      collapse = ""), output = "character")
    
    label_error <- paste(c("RMSE = ",tibble_error$RMSE_mean, "\n",
                           "RRMSE = " , 100*tibble_error$RRMSE_mean, "% \n",
                           "EF = ", 100*tibble_error$EF_mean, "% \n",
                           "%Var imp = ", round(variance_decomposition(data_fitted_observed)$prop_var_imputation, digits =2 )*100, "%"),
                         collapse = "")
    
    plot_output <- plot_output + 
      annotate(
      x = min(data_fitted_observed$observed),
      y = max(data_fitted_observed$fitted),
      geom = "label",
      label = label_error,
      hjust = 0,
      vjust = 0.8, 
      # parse = TRUE,
      fill = as.character(theme_get()$panel.background$fill),
      size = (size_text - 6)/ggplot2::.pt  # transform pt in mm
    ) 
  }
  
  return(plot_output)
  
}



my_scale_x <- function(){
  
  return(scale_x_discrete(labels=c('diff_slope_biomass'=parse(text = TeX('$\\Delta \\mu_{biom}$')),
                                   "diff_IC_SC_asymp_height_C" = parse(text = TeX("$\\Delta_{IC-SC, C, height}$")),
                                   "diff_IC_SC_asymp_height_L" = parse(text = TeX("$\\Delta_{IC-SC, L, height}$")),
                                   "diff_IC_SC_max_LAI_C" = parse(text = TeX("$\\Delta_{IC-SC, C, LAI}$")),
                                   "diff_IC_SC_max_sla_GLT_C" = parse(text = TeX("$\\Delta_{IC-SC, C, SLA}$")),
                                   "diff_IC_SC_max_sla_GLT_L" = parse(text = TeX("$\\Delta_{IC-SC, L, SLA}$")),
                                   "diff_IC_SC_max_LAI_L" = parse(text = TeX("$\\Delta_{IC-SC, L, LAI}$")),
                                   'diff_lambda_height'=parse(text = TeX('$\\Delta_{\\lambda}  height$')),
                                   'diff_lambda_biomass'=parse(text = TeX('$\\Delta_{\\lambda}  biom$')),
                                   'diff_max_LAI' = parse(text = TeX('$\\Delta_{max, LAI}$')), 
                                   'diff_slope_height'=parse(text = TeX('$\\Delta_{\\mu, height}$')),
                                   'diff_asymp_height'=parse(text = TeX('$\\Delta_{max, height}$')),
                                   "diff_IC_SC_slope_height_L" = parse(text = TeX("$\\Delta_{IC-SC, \\mu, height, L}$")),
                                   "diff_IC_SC_slope_height_C" = parse(text = TeX("$\\Delta_{IC-SC} \\mu, height, C$")),
                                   "diff_IC_SC_slope_biomass_L" = parse(text = TeX("$\\Delta_{IC-SC,\\mu, biom, L}$")),
                                   "diff_IC_SC_slope_biomass_C" = parse(text = TeX("$\\Delta_{IC-SC, biom, \\mu,  C}$")),
                                   "diff_IC_SC_lambda_height_L" = parse(text = TeX("$\\Delta_{IC-SC, L, \\lambda, height}$")),
                                   "diff_IC_SC_lambda_height_C" = parse(text = TeX("$\\Delta_{IC-SC, C, \\lambda, height}$")),
                                   "diff_IC_SC_lambda_biomass_L" = parse(text = TeX("$\\Delta_{IC-SC, L, \\lambda, biom}$")),
                                   "diff_IC_SC_lambda_biomass_C" = parse(text = TeX("$\\Delta_{IC-SC, C, \\lambda, biom}$")),
                                   'cover_before_800_normalized_tt' =  'cover',
                                   'cover_before_800_normalized_tt_SC_C' = parse(text = TeX("$\\Delta_{IC-SC, C, cover}$")),   
                                   'cover_before_800_normalized_tt_SC_L' = parse(text = TeX("$\\Delta_{IC-SC, L, cover}$")),      
                                   "density_factor" = "density")))
  
  
}



my_scale_fill_correlations <- function( col_correlations =c("#EFA86E","#db735c", "#555555")){
  return(scale_fill_manual(values = c("POS"  = col_correlations[1],
                                      "NEG" = col_correlations[2],
                                      "NA"=col_correlations[3]),
                           limits = force))
}




# Plot y vs numeric covariates
function_ggplot_numeric <- function(data, response = yield_cereal, n_breaks = 4){
  
  
  vec_y = data %>% distinct(response) %>% pull
  
  breaks_y = round(seq(from = min(vec_y), to = max(vec_y), length.out = n_breaks), digits = 1)
  
  plot_output <- data %>%
    distinct %>% 
    ggplot(aes(x = covariate, y=  {{response}})) +
    geom_point(alpha = .4) + 
    geom_line(stat = 'smooth', alpha = .3)+
    geom_ribbon(stat='smooth', se=TRUE, alpha=0.1,  #SO to get nice confidence intervals  https://stackoverflow.com/questions/29235114/control-transparency-of-smoother-and-confidence-interval
                aes(color = NULL)) + 
    scale_y_continuous( breaks = breaks_y) + 
    theme(axis.title = element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank(), plot.margin =  unit(c(0,0,0,0), "cm") )
  
  return(plot_output)
}

# Plot y vs quali covariates
function_ggplot_qualitative <- function(data, response = yield_cereal,  width_jitter_points = 0.05, n_breaks = 4){
  
  vec_y = data %>% distinct(response) %>% pull
  
  breaks_y = round(seq(from = min(vec_y), to = max(vec_y), length.out = n_breaks), digits = 1)
  
  plot_output <- data %>%
    distinct %>% 
    mutate(covariate = forcats::fct_reorder(covariate, {{response}}, .desc = TRUE)) %>%
    ggplot(aes(x = covariate, y=  {{response}})) +
    geom_point(alpha = .3, position = position_jitter(width = width_jitter_points))   +
    geom_boxplot(alpha = .2)+
    scale_y_continuous(breaks = breaks_y) +
    ylim(c(0, max(data %>% select({{response}}) %>% pull) + 1)) +
    # theme(axis.title = element_blank(), axis.ticks = element_blank(),
    # axis.text.x = element_text(vjust = 1,
    # angle = 45,
    # margin = ggplot2::margin(t = -12, b = 12)), plot.margin =  unit(c(0,0,0,0), "cm"))
    theme(axis.title = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), plot.margin =  unit(c(0,0,0,0), "cm"))
    
  
  return(plot_output)
}




# Plot the importance and add marginal plots witin them
plot_importance <- function(data_models, threshold = 10,
                            lab_y =TeX("Importance ($\\Delta$ MSE after shuffling covariate)"), size_text = 15,
                            add_marginal_plots = TRUE, y,
                            n_breaks = 4) {
  
  
  data_importance <- get_importance(data_models)
  
  
  data_plot <- data_importance %>% 
    filter(n >= threshold) 
  
  
  
  plot_output <- data_plot %>%
    ggplot(aes(x=variable, mean_importance)) +
    geom_bar(stat = "identity", alpha = .4, aes(fill = sign_cor)) +
    geom_errorbar(aes(ymin = mean_importance -  sd_importance, ymax = mean_importance + sd_importance), width = .5 ) + 
    ylab(lab_y) +
    xlab(NULL) +
    labs(fill = TeX('Sign $\\tau$ Kendall'), title = NULL) +
    my_scale_fill_correlations() +
    my_scale_x()+
    theme(legend.position = "bottom",
          axis.text.y = element_text(angle = 60, hjust = 0.5),
          title = element_text(size = size_text),
          legend.text = element_text(size = size_text - 2 ),
          legend.title = element_text(size = size_text - 2 ),
          axis.title = element_text(size = size_text - 2),
          axis.text = element_text(size = size_text +2)) +
    coord_flip() 
  
  
  
  
  if(add_marginal_plots){
    
    
    marginal_plots_quali <- data_models %>%
      unnest(data) %>%
      select(-c(n_iterations, mod.gmerf, .imp, experiment_id, management) & where(is.character), response = any_of(y))%>%
      pivot_longer(-  response , names_to = "variable", values_to = "covariate")  %>%
      mutate(variable = str_replace(variable, "cereal" , "C")) %>%  
      mutate(variable = str_replace(variable, "legume" , "L")) %>% 
      group_by(variable) %>% 
      nest() %>% 
      mutate(marginal_plot = map(data, function_ggplot_qualitative, response = response, n_breaks = n_breaks)) %>% 
      ungroup
    
    
    
    marginal_plots_numeric <- data_models %>%
      unnest(data) %>%
      select(-c(n_iterations, mod.gmerf, .imp, experiment_id, management) & where(is.numeric), response =any_of(y))%>%
      select(-any_of(y)) %>% 
      pivot_longer(-  response , names_to = "variable", values_to = "covariate") %>% 
      mutate(variable = str_replace(variable, "cereal" , "C")) %>%  
      mutate(variable = str_replace(variable, "legume" , "L")) %>%  
      group_by(variable) %>% 
      nest() %>% 
      mutate(marginal_plot = map(data, function_ggplot_numeric, response = response, n_breaks = n_breaks )) %>% 
      ungroup
    
    
    
    data_marginal_plots <- bind_rows(marginal_plots_numeric, marginal_plots_quali) %>% 
      mutate(variable = str_replace(variable, "cereal" , "C")) %>%  
      mutate(variable = str_replace(variable, "legume" , "L")) %>%  
      inner_join(data_importance ) %>%
      filter(n >= threshold) %>%
      select(-data)%>% 
      mutate(variable = fct_reorder(variable, mean_importance)) %>% 
      arrange(desc(variable))
    
    plot_inset <- cowplot::plot_grid(plotlist = data_marginal_plots$marginal_plot, ncol =1)
    
    
    width_insets_plot <- 1.2* max(data_marginal_plots$mean_importance + data_marginal_plots$sd_importance)
    # max_y <- max(data_marginal_plots$mean_importance + data_marginal_plots$sd_importance )+ 0.5/max(data_marginal_plots$mean_importance) + width_insets_plot
    max_y <- max(data_marginal_plots$mean_importance + data_marginal_plots$sd_importance )+ 0.25*max(data_marginal_plots$mean_importance) + width_insets_plot
    
    plot_output <- plot_output+
      ylim(c(0, max_y)) +
      annotation_custom(grob = ggplotGrob(plot_inset), ymax = max_y, ymin = max_y - width_insets_plot)
    
    
  }
  
  return(plot_output)
  
}



