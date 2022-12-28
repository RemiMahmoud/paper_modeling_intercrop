
#' evaluate efficiency criterion of fitted vs simulated dataset 


evaluate_error <- function (data, observed = "observed", fitted = "fitted") {
  metrics <- data %>%
    select(observed =  all_of(observed), fitted = all_of(fitted)) %>% 
    drop_na() %>%
    summarise(n = n(),
              mean_observed = mean(observed), 
              mean_fitted = mean(fitted),
              bias = mean(observed - fitted),
              bias_squared = bias^2,
              SSE = sum((observed - 
                           fitted)^2),
              MSE = SSE/n,
              RMSE = MSE^0.5, 
              RRMSE = RMSE/mean_observed, 
              RRMSE_IQR = RMSE/IQR(observed),
              RRMSE_SD = RMSE/sd(observed),
              MAE = mean(abs(observed - fitted)),
              RMAE = MAE/mean(abs(observed)),
              RMAE_IQR = MAE/IQR(observed),
              RMAEP = mean(abs(observed - fitted)/abs(observed)), 
              EF = 1 - sum((observed - fitted)^2)/sum((observed - 
                                                         mean_observed)^2),
              index_willmott = 1 - sum((observed - fitted)^2)/sum((abs(fitted - mean(observed)) + 
                                                                     abs(observed - mean(observed)))^2), 
              SDSD = (sd(fitted) - sd(observed))^2 * (n - 1)/n,
              LCS = 2 * sd(observed) * sd(fitted) * (1 - cor(observed, fitted)) * 
                (n - 1)/n,
              NU = (1 - (cov(observed, fitted)/var(fitted)))^2 * 
                var(fitted) * (n - 1)/n, 
              LC = (1 - cor(observed, 
                            fitted)^2) * var(observed) * (n - 1)/n, 
              r_pearson = cor(fitted, 
                              observed, method = "pearson"),
              p_pearson = ifelse(n > 
                                   2, cor.test(fitted, observed, method = "pearson")$p.value, 
                                 NA), 
              r_kendall = cor(fitted, observed, method = "kendall"), 
              p_kendall = ifelse(n > 2, cor.test(fitted, observed, 
                                                 method = "kendall")$p.value, NA),
              r_squared = cor(fitted, 
                              observed, method = "pearson")^2)
  return(metrics)
}
