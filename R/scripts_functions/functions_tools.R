
#' function to split data in a same way according to different models


make_data_split <- function(data, frac_train = 0.75, set_seed = TRUE){
  
  if(set_seed){ set.seed(as.integer(Sys.Date()))} # Set a different seed for each day
  # data_to_split <- data 
  data_splitted <- rsample::initial_split(data, prop = frac_train)
  
  return(data_splitted)
}

#' function to add quickliy y =x line
geom_y_x <- function(){geom_abline(slope =1, intercept = 0) } 


custom_coerce_fct_na_explicit <- function(x){
  if (is.factor(x) & anyNA(x)) {
    forcats::fct_explicit_na(x, na_level = "NA")
  } else {
    x
  }
}

custom_gg_miss_fct <- function(x, fct){
  
  fct <- rlang::enquo(fct)
  
  data <- x %>%
    # protect against error where grouping by missing value leads to
    # warning message from dplyr about explicit
    dplyr::mutate_at(vars(!!fct), .funs = custom_coerce_fct_na_explicit) %>%
    dplyr::group_by(!!fct) %>%
    miss_var_summary()
  
  ggobject <-
    ggplot(data,
           aes_string(x = quo_name(fct),
                      y = "variable",
                      fill = "pct_miss")) +
    geom_tile() +
    viridis::scale_fill_viridis(name = "% Miss") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  
  return(ggobject)
}


shorten_exp_id <- function(data){
  return(data%>%
           mutate(experiment_id = str_replace_all(experiment_id, "Sanmarco", "sanmarco")) %>%
           mutate(experiment_id = str_replace_all(experiment_id, "sanmarco", "SM")) %>%
           mutate(experiment_id = str_replace_all(experiment_id, "SMargentano_SMargentano", "SMA")) %>%
           mutate(experiment_id = str_replace_all(experiment_id, "devad", "")) %>%
           mutate(experiment_id = str_replace_all(experiment_id, "Angers", "Ang")) %>%
           mutate(experiment_id = str_replace_all(experiment_id, "Auzeville", "Auz")) %>%
           mutate(experiment_id = str_replace_all(experiment_id, "Copenhagen", "Copen")) %>% 
           mutate(experiment_id = str_replace_all(experiment_id, "Reading_reading", "Reading")) %>%
           mutate(experiment_id = str_replace_all(experiment_id, "_taastrup", "")) %>% 
           mutate(experiment_id = str_replace_all(experiment_id, "_jyn", "")) %>% 
           mutate(experiment_id = str_replace_all(experiment_id, "_hbg", "")) %>% 
           mutate(experiment_id = str_replace_all(experiment_id, "_TO", "")) %>% 
           mutate(experiment_id = str_replace_all(experiment_id, "_tesgues", "")) %>% 
           mutate(experiment_id = str_replace_all(experiment_id, "_inra", "")) %>% 
           mutate(experiment_id = str_replace_all(experiment_id, "_les_roches", ""))  %>%
           mutate(experiment_id = str_replace_all(experiment_id, "marinette", "marin"))  %>%
           mutate(experiment_id = str_replace_all(experiment_id, "Kassel_kassel", "Kass"))  %>%
           mutate(experiment_id = str_replace_all(experiment_id, "thorigne", "thor")))
}
