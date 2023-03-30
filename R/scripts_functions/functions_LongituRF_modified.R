
# Taken from https://github.com/sistm/LongituRF/blob/master/R/LongituRF.R

MERF_modified <- function(X,Y,id,Z,iter=100,ntree=500, time, sto = "none", delta = 0.001, verbose = TRUE){
  q <- dim(Z)[2]
  nind <- length(unique(id))
  btilde <- matrix(0,nind,q) #### Pour la ligne i, on a les effets al?atoires de l'individu i
  sigmahat <- 1 #### init
  Btilde <- diag(rep(1,q)) ### init
  epsilonhat <- rep(0,length(Y))
  id_btilde <- unique(id)
  Tiime <- sort(unique(time))
  omega <- rep(0,length(Y))
  sigma2 <- 1
  Vrai <- NULL
  inc <- 1
  OOB <- NULL
  
  
  ########## ADDENDUM REMI MAHMOUD ###############################################################
  #########################################################################################################
  #########################################################################################################
  
  xnam=names(X)
  
  features_selected <- xnam # All features selected at first iteration
  
  #########################################################################################################
  #########################################################################################################
  
  
  if (class(sto)=="character"){
    
    if ( sto=="none"){
      for (i in 1:iter){
        ystar <- rep(NA,length(Y))
        for (k in 1:nind){ #### on retrace les effets al?atoires
          indiv <- which(id==unique(id)[k])
          ystar[indiv] <- Y[indiv]- Z[indiv,,drop=FALSE]%*%btilde[k,]
        }
        
        
        
        
        ########## ADDENDUM REMI MAHMOUD ###############################################################
        #########################################################################################################
        #########################################################################################################
        
        forest.formula=as.formula(paste("ystar ~ ", paste(features_selected, collapse= "+")))
        
        boruta.data=cbind(ystar, X[features_selected])
        mod.boruta <- Boruta(forest.formula, boruta.data)
        
        # RF FIT
        features_selected <- getSelectedAttributes(mod.boruta)
        
        forest <- randomForest(X[features_selected],ystar,ntree=ntree, importance = TRUE) ### on construit l'arbre
        
        #########################################################################################################
        #########################################################################################################
        
        # forest <- randomForest(X,ystar,mtry=mtry,ntree=ntree, importance = TRUE) ### on construit l'arbre
        fhat <- predict(forest)
        OOB[i] <- forest$mse[ntree]
        for (k in 1:nind){
          indiv <- which(id==unique(id)[k])
          V <- Z[indiv,, drop=FALSE]%*%Btilde%*%t(Z[indiv,, drop=FALSE])+diag(as.numeric(sigmahat),length(indiv),length(indiv))
          btilde[k,] <- Btilde%*%t(Z[indiv,, drop=FALSE])%*%solve(V)%*%(Y[indiv]-fhat[indiv])
          epsilonhat[indiv] <- Y[indiv] -fhat[indiv] -Z[indiv,, drop=FALSE]%*%btilde[k,]
        }
        
        sigm <- sigmahat
        sigmahat <- sig(sigma = sigmahat,id = id, Z = Z, epsilon = epsilonhat,Btilde = Btilde)
        Btilde  <- bay(bhat = btilde,Bhat = Btilde,Z = Z,id = id,sigmahat = sigm)
        Vrai <- c(Vrai, logV(Y,fhat,Z,time,id,Btilde,0,sigmahat,sto))
        if (i>1) inc <- abs((Vrai[i-1]-Vrai[i])/Vrai[i-1])
        
        
        
        
        if(verbose) {
          if(i%%10 ==0) {
            print(paste0("Iteration number ", i) )
            print(paste0("Conv criterion: ", round(inc, digits = 3)))}
        }
        
        if (inc < delta) {
          print(paste0("stopped after ", i, " iterations."))
          sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance = Vrai,id=id, time=time, OOB =OOB)
          class(sortie) <- "longituRF"
          return(sortie)
        }
      }
      sortie <- list(forest=forest,random_effects=btilde,var_random_effects=Btilde,sigma=sigmahat, id_btilde=unique(id), sto= sto, vraisemblance=Vrai,id=id, time=time, OOB =OOB)
      class(sortie) <- "longituRF"
      return(sortie)
    }
  }
}


predict.longituRF <- function(object, X,Z,id,time,...){
  n <- length(unique(id))
  id_btilde <- object$id_btilde
  f <- predict(object$forest,X)
  Time <- object$time
  id_btilde <- object$id_btilde
  Ypred <- rep(0,length(id))
  id.app=object$id
  if (object$sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      k <- which(id_btilde==unique(id)[i])
      Ypred[w] <- f[w] + Z[w,, drop=FALSE]%*%object$random_effects[k,]
    }
    return(Ypred)
  }
  return(Ypred)
}


#' Title
#'
#'
#' @import stats
#'
#' @keywords internal
sig <- function(sigma,id,Z, epsilon, Btilde){ #### fonction d'actualisation du param?tre de la variance des erreurs
  nind <- length(unique(id))
  Nombre <- length(id)
  sigm <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
    sigm <- sigm + t(epsilon[w])%*%epsilon[w] + sigma*(length(w)-sigma*(sum(diag(solve(V)))))
  }
  sigm <- sigm/Nombre
  return(sigm)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
bay <- function(bhat,Bhat,Z,id, sigmahat){ #### actualisation des param?tres de B
  nind <- length(unique(id))
  q <- dim(Z)[2]
  Nombre <- length(id)
  D <- 0
  for (j in 1:nind){
    w <- which(id==unique(id)[j])
    V <- Z[w,, drop=FALSE]%*%Bhat%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
    D <- D+ (bhat[j,]%*%t(bhat[j,]))+ (Bhat- Bhat%*%t(Z[w,, drop=FALSE])%*%solve(V)%*%Z[w,, drop=FALSE]%*%Bhat)
  }
  D <- D/nind
  return(D)
}

#' Title
#'
#' @import stats
#'
#' @keywords internal
# Moy <- function(id,Btilde,sigmahat,Phi,Y,Z){
#   S1<- 0
#   S2<- 0
#   nind <- length(unique(id))
#   for (i in 1:nind){
#     w <- which(id==unique(id)[i])
#     V <- Z[w,, drop=FALSE]%*%Btilde%*%t(Z[w,, drop=FALSE])+diag(as.numeric(sigmahat),length(w),length(w))
#     S1 <- S1 + t(Phi[w,, drop=FALSE])%*%solve(V)%*%Phi[w,, drop=FALSE]
#     S2 <- S2 + t(Phi[w,, drop = FALSE])%*%solve(V)%*%Y[w]
#   }
#   return(solve(S1)%*%S2)
# }


logV <- function(Y,f,Z,time,id,B,gamma,sigma, sto){
  Vraisem <- 0
  if (sto=="none"){
    for (i in 1:length(unique(id))){
      w <- which(id==unique(id)[i])
      V <- Z[w,,drop=FALSE]%*%B%*%t(Z[w,,drop=FALSE])+diag(as.numeric(sigma),length(w),length(w))
      V2 = log(det(V))+ t(Y[w]-f[w])%*%solve(V)%*%(Y[w]-f[w])
      if (V2<Inf){
        Vraisem <- Vraisem + V2
      }
    }
    return(Vraisem)
  }
}
