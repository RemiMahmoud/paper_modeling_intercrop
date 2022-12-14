#' Mixed Random Forest
#'
#' The function to fit a random forest with random effects.
#'
#' param Y The outcome variable.
#' param X A data frame or matrix contains the predictors.
#' param random A string in lme4 format indicates the random effect model.
#' param data The data set as a data frame.
#' param initialRandomEffects The initial values for random effects.
#' param ErrorTolerance The tolerance for log-likelihood.
#' param MaxIterations The maximum iteration times.
#'
#' return A list contains the random forest ($forest), mixed model ($MixedModel), and random effects ($RandomEffects).
#' See the example below for the usage.

#' export
#' import randomForest lme4
#' examples
#'
#' # load the sleepstudy data set from the lme4 package
#' library(lme4)
#' data(sleepstudy)
#'
#' tmp = MixRF(Y=sleepstudy$Reaction, X=as.data.frame(sleepstudy$Days), random='(Days|Subject)',
#'             data=sleepstudy, initialRandomEffects=0, ErrorTolerance=0.01, MaxIterations=100)
#'
#' # tmp$forest
#' 
#' # tmp$MixedModel
#' 
#' # tmp$RandomEffects


MixRF <- function(Y, X, random, data, initialRandomEffects = 0, ErrorTolerance = 0.001, MaxIterations = 1000, 
                 importance = FALSE, ntree = 500, mtry = max(floor(ncol(X)/3), 1), nodesize = 5, maxnodes = NULL, verbose= TRUE) {
  
  Target = Y
  
  # Condition that indicates the loop has not converged or run out of iterations
  ContinueCondition = TRUE
  
  iterations <- 0
  
  # Get initial values
  AdjustedTarget <- Target - initialRandomEffects
  oldLogLik <- -Inf
  
  
  conv_crit <- rep(0, MaxIterations)
  data_b <- as.data.frame(matrix(rep(0,MaxIterations*length(unique(data$experiment_id))), ncol = length(unique(data$experiment_id))))
  names(data_b) <- unique(data$experiment_id)
  
  xnam <- names(X)
  features_selected <- xnam # All features selected at first iteration
  
  while(ContinueCondition){
    
    iterations <- iterations+1
    
    # randomForest
    
    forest.formula=as.formula(paste("AdjustedTarget ~ ", paste(features_selected, collapse= "+")))
    
    # FEATURE SELECTION RANDOM FOREST
    boruta.data=cbind(AdjustedTarget, X[features_selected])
    mod.boruta <- Boruta(forest.formula, boruta.data)
    
    # RF FIT
    features_selected <- getSelectedAttributes(mod.boruta)
    
    rf = randomForest(X[features_selected], AdjustedTarget, 
                      importance = importance, ntree = ntree, mtry = mtry, nodesize = nodesize, maxnodes = maxnodes)
    
    # y - X*beta (out-of-bag prediction)
    resi = Target - rf$predicted
    
    ## Estimate New Random Effects and Errors using lmer
    f0 = as.formula(paste0('resi ~ -1 + ',random))
    lmefit <- lmer(f0, data=data)
    
    
    # data_b[iterations,] = bi
    
    # check convergence
    newLogLik <- as.numeric(logLik(lmefit))
    
    ContinueCondition <- (abs(newLogLik-oldLogLik)>ErrorTolerance & iterations < MaxIterations)
    
    
    
    conv_crit[iterations] <- abs(newLogLik-oldLogLik)
    
    
    if(verbose) {
      if(iterations%%10 ==0) {
        print(paste0("Iteration number ", iterations, "\n") )
        print(paste0("Conv criterion: ",round(abs(newLogLik-oldLogLik), digits = 4)))}
    }
    
    oldLogLik <- newLogLik
    
    # Extract random effects to make the new adjusted target
    AllEffects <- predict(lmefit)
    
    #  y-Zb
    AdjustedTarget <- Target - AllEffects
    
    
    
  }
  
  result <- list(forest=rf, MixedModel=lmefit, RandomEffects=ranef(lmefit),
                 IterationsUsed=iterations, conv_crit = conv_crit)
  
  return(result)
}
