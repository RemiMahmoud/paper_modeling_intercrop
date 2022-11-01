
#' evaluate efficiency criterion of fitted vs simulated dataset 

efficiency <- function(data, observed = "observed", simulated = "simulated") 
{
  EF <- data %>%
    select(observed = observed, simulated = simulated) %>% 
    drop_na() %>%
    summarise( mean_observed = mean(observed), 
               EF = 1 - sum((observed - simulated)^2)/sum((observed - 
                                                             mean_observed)^2)) %>% 
    pull(EF)
  return(EF)
  
}