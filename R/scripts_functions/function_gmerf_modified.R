
gmerf_modified= function (y, cov, group, xnam=NULL, znam=NULL, family='binomial', bizero=NULL,
                          itmax=30, toll=0.02, weights_lme = NULL, control_lme = list(), verbose = TRUE ) {
  #argomenti:
  #-y=vettore con le risposte
  #-cov=data frame con le covariate di ogni unit� statistica
  #-group=vettore factor che, per ogni unit� statistica, dice il gruppo a cui appartiene
  #-xnam=vettore coi nomi delle covariate da usare nella random forest
  #-znam=vettore coi nomi delle covariate da usare nei random effects
  #-bizero=matrice in cui ogni colonna contiene i coefficienti dei random effects
  #	   il primo valore di ogni colonna � l'intercept, gli altri le covariate znam
  
  # verbose= TRUE se dobbiamo mostraro il numero di iterazioni
  
  #assumo che group e bizero siano coerenti, cio� b[,i] corrisponda a levels(group)[i]
  
  
  ######################################
  ####	STEP 1: Inizializzazione  ######
  ######################################
  
  N <- length(y) #numero di osservazioni
  n=length(levels(group)) #numero di gruppi
  q <- length(znam)+1	# numero di covariate + random intercept
  Zi=NULL
  z.not.null=!(is.null(znam)) #controllo se ci sono covariate incluse nei random effects
  if (z.not.null) Zi=cov[znam] #covariate dei random effects
  Zi.int= cbind(rep(1,N),Zi) #random intercept + covariate dei random effects
  
  #Inizializzo (se � NULL) bi a 0
  if( is.null(bizero) ){
    bi <- NULL
    for(i in 1:n) bi=cbind(bi,rep(0,q))
  }
  if( !is.null(bizero) ) bi=bizero
  lev=levels(group) #nomi dei gruppi
  bi=data.frame(bi)
  names(bi)=lev
  all.bi=list()  #i b_i di ogni iterazione
  all.bi[[1]]=bi
  
  #se xnam � NULL assumo che tutte le variabili di cov siano da usare
  if(is.null(xnam)) xnam=names(cov)
  
  #group deve essere un factor, altrimenti da errore
  if(!is.factor(group)) stop('Argomento "group" deve essere un factor')
  
  if(z.not.null){
    glmer.formula=as.formula(paste("y ~ ( 1+", paste(znam, collapse= "+"), " | group )"))
    lme.random.formula = as.formula(paste(" ~ 1+", paste(znam, collapse= "+"), " | group "))
  }
  if(!z.not.null){
    glmer.formula=as.formula(paste("y ~ ( 1 | group) "))
    lme.random.formula = as.formula(paste(" ~ 1|group "))
  }
  
  
  # forest.formula=as.formula(paste("target ~ ", paste(xnam, collapse= "+")))
  
  
  ##################################################
  ####	STEP 2: GLM per inizializzare mu_ij  #######
  ##################################################
  
  glm.formula=as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  glm.data= cbind(y,cov[xnam])
  glm.0=glm(glm.formula, data=glm.data, family=family)
  mu.ij.0=glm.0$fitted.values
  eta.est=glm.0$family$linkfun(mu.ij.0) #stime del glm
  linkf=glm.0$family$linkfun
  linkinv=glm.0$family$linkinv
  
  
  ####################################################
  ####	STEP 3-4: Stima iterativa del modello  #######
  ####################################################
  
  
  it=1
  converged=FALSE
  
  features_selected <- xnam # All features selected at first iteration
  forest.formula=as.formula(paste("target ~ ", paste(features_selected, collapse= "+")))
  
  
  conv_crit <- rep(0, itmax)
  data_b <- as.data.frame(matrix(rep(0,itmax*length(lev)), ncol = length(lev)))
  names(data_b) <- lev
  
  while(!converged && it<itmax) {
    
    
    #random forest
    target=rep(0,N) #target=eta-Z%*%b
    for (i in 1:N) {
      b.temp=as.matrix(bi[group[i]], nrow=q, ncol=1)
      z.temp=as.matrix(Zi.int[i,], nrow=1, ncol=q)
      target[i]= eta.est[i] - z.temp%*%b.temp
    }
    
    
    
    # FEATURE SELECTION RANDOM FOREST
    boruta.data=cbind(target, cov[features_selected])
    mod.boruta <- Boruta(forest.formula, boruta.data)
    
    # RF FIT
    features_selected <- getSelectedAttributes(mod.boruta)
    forest.data=cbind(target, cov[features_selected])
    
    # if(length(features_selected) <2){print(features_selected)}
    forest.formula=as.formula(paste("target ~ ", paste(features_selected, collapse= "+")))
    # if(length(features_selected) <2){print(forest.formula)}
    
    
    forest=randomForest(forest.formula, forest.data)
    
    f.x_ij=forest$predicted
    
    #glm con mixed effects
    #glmer.data=cbind(y,group)
    #if(z.not.null) glmer.data=cbind(glmer.data, Zi)
    #glmer.data=data.frame(glmer.data) #altrimenti glmer non funziona
    glmer.data = data.frame(y,group)
    lme.data = data.frame(y = y - f.x_ij, group) # include offset in the data frame for lme
    
    if(z.not.null) glmer.data=data.frame(glmer.data, Zi)
    
    if(family == "gaussian"){lme.fit = lme(data = lme.data, random = lme.random.formula, fixed = as.formula("y~1"), weights = weights_lme, control = control_lme)}
    else
    {glmer.fit= glmer(glmer.formula, glmer.data, family=family, offset=f.x_ij)}
    
    
    
    #voglio mantenere l'ordine degli elementi dei b_i
    select=c("(Intercept)",znam) #tutti i b_i da estrarre
    
    if(family == "gaussian"){lme.bi =ranef(lme.fit)[select]}else{glmer.bi=ranef(glmer.fit)$group[select]}
    
    
    
    #convergenza dei b_i
    bi.old=bi
    if(family == "gaussian"){bi=data.frame(t(lme.bi))} else { bi=data.frame(t(glmer.bi))}
    
    
    
    names(bi)=lev
    data_b[it,] = bi
    
    # if(verbose){print(bi)}
    
    
    diff.t=abs(bi.old-bi)
    n.diff=max(diff.t) #uso la norma infinito(max)
    ind=which(diff.t==n.diff, arr.ind=T)
    n.old=abs(bi.old[ind])
    converged= n.diff/n.old <toll
    
    
    if(verbose) {
      if(it%%10 ==0) {
      print(paste0("Iteration number ", it, "\n") )
      print(paste0("Conv criterion: ", n.diff/n.old))}
    }
    
    
    conv_crit[it] = n.diff/n.old
    it=it+1
    all.bi[[it]]=bi
  }
  
  ###############################################
  ####	STEP 5: Preparazione output  ############
  ###############################################
  
  #se non ho convergenza do un messaggio di errore
  if(!converged) {
    warning('Numero massimo di iterazioni superato, non si � arrivati a convergenza')
  }
  
  if(family == "gaussian"){mixed.model = lme.fit} else {mixed.model = glmer.fit}
  result=list(mixed.model,forest,bi,it,converged,all.bi,linkf,linkinv,xnam,znam,family, f.x_ij, features_selected, conv_crit, data_b[-nrow(data_b),])
  names(result)=c('mixed.model', 'forest.model', 'rand.coef', 'n.iteration',
                  'converged','all.rand.coef','linkf','linkinv','forest.var',
                  'random.eff.var','family', 'f.x_ij', 'features_selected_boruta', "convergence_criterion", "data_b")
  class(result)='gmerf'
  result
}






fitted.gmerf=function(gm, type='response', alpha=0.5) {
  #i fitted values corretti sono gi� quelli della funzione glmer, poich�
  #ho incorporato i fitted values della random forest nel modello
  mu=fitted(gm$mixed.model)
  #errore se il tipo di risposta non � tra questi
  allowed=c('response', 'mu', 'eta')
  msg=paste('Possible choices are ', allowed[1],', ' ,allowed[2],' and ',
            allowed[3], sep='')
  if(sum(type==allowed)==0) stop('Type of prediction not available:',msg)
  if(type=='mu') {ans = mu + gm$f.x_ij}
  if(gm$family=='binomial' && type=='response') ans=mu>alpha
  if(type=='eta') ans=gm$linkf(mu)
  ans
}


summary.gmerf=function(gm) {
  print('Mixed effects model') #summary del mixxed effects model
  print(summary(gm$mixed.model))
  str=ifelse(gm$converged, 'Converged', 'Did not converge')
  print(paste(str , 'after', gm$n.iteration, 'iterations')) #dice se c'� convergenza
}





predict.gmerf=function(gm, newdata, group, type='response', alpha=0.5,
				predict.all=FALSE, re.form=NULL, newparam=NULL,
				terms=NULL, allow.new.levels=TRUE, na.action=na.pass,
				random.only=FALSE) {

#chiamo semplicemente i metodi predict di gmerf e randomForest
	forest.data=newdata[gm$forest.var]
	glmer.data=data.frame(newdata[gm$random.eff.var],group)
	p1=predict(gm$glmer.model,newdata=glmer.data,newparam=newparam,re.form=re.form,
			random.only=random.only, type='link', na.action=na.action,
			allow.new.levels=allow.new.levels)
	p2=predict(gm$forest.model,forest.data,predict.all=predict.all)

	#errore se il tipo di risposta non � tra questi
	allowed=c('response', 'mu', 'eta')
	msg=paste('Possible choices are ', allowed[1],', ' ,allowed[2],' and ',
			allowed[3], sep='')
	if(sum(type==allowed)==0) stop('Type of prediction not available:',msg)

	eta=p1+p2
	mu=gm$linkinv(eta)
	if(type=='eta') ans=eta #eta:predico eta
	if(type=='mu') ans=mu #predico la probabilit�
	if(gm$family=='binomial' && type=='response') ans=mu>alpha #predico la risposta basandomi su alpha
	ans
}
