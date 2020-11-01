IEHC = function(data,X,G,M,
                chisq_app = "3M",
                combinations_return=TRUE,
                combination_preference = "All", # other options: "optim", "adapt", "Fisher","Burden","KM","ACAT"
                acc = 5e-10,
                accurate_app_threshold = -log10(0.05), # if the min_pvalue<0.05, the accurate approximation on min_quantile in optim combination is used instead of approximation with Liu's method.
                max_core = 4
){

library(survival)
library(parallel)
ACAT<-function(Pvals,Weights=NULL){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
        stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
        stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(sum(Pvals==0)>=1)
    is.one<-(sum(Pvals==1)>=1)
    if (is.zero && is.one){
        stop("Cannot have both 0 and 1 p-values!")
    }
    if (is.zero){
        return(0)
    }
    if (is.one){
        warning("There are p-values that are exactly 1!")
        return(1)
    }

    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(Weights)){
        Weights<-rep(1/length(Pvals),length(Pvals))
    }else if (length(Weights)!=length(Pvals)){
        stop("The length of weights should be the same as that of the p-values")
    }else if (sum(Weights<0)>0){
        stop("All the weights must be positive!")
    }else{
        Weights<-Weights/sum(Weights)
    }


    #### check if there are very small non-zero p values
    is.small<-(Pvals<1e-16)
    if (sum(is.small)==0){
        cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
    }else{
        cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
        cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
    }
    #### check if the test statistic is very large.
    if (cct.stat>1e+15){
        pval<-(1/cct.stat)/pi
    }else{
        pval<-1-pcauchy(cct.stat)
    }
    return(pval)
}

WeightSumChisq = function(weight,df=rep(1,length(weight))){
  # Calculate the df and the noncentral parameter of an approximated noncentral chisq distribution for the weighted sum of independent chisq distributions.
  # Option weight: a vector of weights for each chisq.
  # Option df: a vector of df for each chisq.
  # Return the df(l) and the noncentrality (\delta) of the approximation noncentral chisq.
  # Reference: Liu (2009). Match the skewness.
  
  # Sum of weights to different powers and terms related to criterions
  S1 = sum(weight)    # equal to c1[1]
  S2 = sum(weight^2)  # equal to c1[2]
  S3 = sum(weight^3)
  S4 = sum(weight^4)  
  tmp_1 = S4/(S2^2)
  tmp_2 = (S3^2)/(S2^3)
  
  # Calculate a=sqrt(l+2\delta), \delta, and l
  if(tmp_1<tmp_2){
    a = 1/(S3/(S2^(3/2))-sqrt(tmp_2-tmp_1))
    delta = a^3*S3/(S2^(3/2))-a^2
    l = a^2-2*delta
  }
  else{
    a = (S2^(3/2))/S3    
    delta = 0
    l = a^2
  }
  #print(c(tmp_1,tmp_2,a,delta,l))
  return(c(mean.Q=S1,SD.Q=sqrt(2*S2),df=l,noncentral=delta,mean.X=l+delta,SD.X=sqrt(2*(l+2*delta))))
}

WeightSumChisq_mod = function(weight,df=rep(1,length(weight))){
  # Calculate the df and the noncentral parameter of an approximated noncentral chisq distribution for the weighted sum of independent chisq distributions.
  # Option weight: a vector of weights for each chisq.
  # Option df: a vector of df for each chisq.
  # Return the df(l) and the noncentrality (\delta) of the approximation noncentral chisq.
  # Reference: Liu (2009). Modified version with matching the kurtosis.
  
  # Sum of weights to different powers and terms related to criterions
  S1 = sum(weight)    
  S2 = sum(weight^2)
  S3 = sum(weight^3)
  S4 = sum(weight^4)  
  tmp_1 = S4/(S2^2)
  tmp_2 = (S3^2)/(S2^3)
  
  # Calculate a=sqrt(l+2\delta), \delta, and l
  if(tmp_1<tmp_2){
    a = 1/(S3/(S2^(3/2))-sqrt(tmp_2-tmp_1))
    delta = a^3*S3/(S2^(3/2))-a^2
    l = a^2-2*delta
  }
  else{
    a = sqrt(1/tmp_1)
    delta = 0
    l = a^2
  }
  #print(c(tmp_1,tmp_2,a,delta,l))
  return(c(mean.Q=S1,SD.Q = sqrt(2*S2),df=l,noncentral=delta,mean.X=l+delta,SD.X=sqrt(2*(l+2*delta))))
}

quan_rho_calculator = function(grid=seq(0,1,by=0.05),df_f,lambda_r,min_pvalue,acc=5e-10,lim=1e4,max_core=4){
  grid_inner = grid[2:(length(grid)-1)]
  numCores = parallel::detectCores()
  
  quan_rho_grid_inner = mclapply(grid_inner,FUN=function(grid_point){
    diff.prob.davies = function(q,tailprob,lambda,df=rep(1,length(lambda)),delta=rep(0,length(lambda)),lim=1e4,acc=5e-10){
      tailprob-CompQuadForm::davies(q,lambda=lambda,lim=lim,acc=acc)$Qq
    }
    quan_tmp = uniroot(diff.prob.davies,lower=0,upper=1e7,tailprob=min_pvalue,lambda=c(grid_point*rep(1,df_f),(1-grid_point)*lambda_r),acc=acc,lim=lim)$root
  },
  mc.cores=min(numCores,max_core))
  quan_rho_grid_inner = as.numeric(unlist(quan_rho_grid_inner))
  names(quan_rho_grid_inner) = as.character(grid_inner)
  return(list(grid=grid,quan_rho_grid_inner=quan_rho_grid_inner))
}

pvalue_minimizePValue_davies = function(u_f,u_r,df_f,lambda_r,min_pvalue=NULL,quan_grid=NULL,acc=5e-10,lim=1e4,ChisqApp="3M",max_core=4){
  if(is.null(min_pvalue)){
    p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
      pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                      lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
      if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,
                                          weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),
                                          method=ChisqApp)$pvalue
      return(pvalue)
    }
    min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)$objective
    min_pvalue = min(c(min_pvalue_inside01,
                       p_rho(0,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc),
                       p_rho(1,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)))
  }
  

  integrand = function(a,min_pvalue,df_f,lambda_r,acc=5e-10,grid=NULL,quan_rho_grid_inner=NULL,max_core=4){
    quan_rho = function(rho,df_f,lambda_r,p_rho_min,u_f){
      if(rho==1) stop("rho is not allowed to be 1.")
      else{
        diff.prob.davies = function(q,tailprob,lambda,df=rep(1,length(lambda)),delta=rep(0,length(lambda)),lim=1e4,acc=5e-10){
          tailprob-CompQuadForm::davies(q,lambda=lambda,lim=lim,acc=acc)$Qq
        }
        quan_tmp = uniroot(diff.prob.davies,lower=0,upper=1e7,tailprob=p_rho_min,lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=lim)$root
        names(quan_tmp) = c()
        quan_tmp_u_r = (quan_tmp - rho*u_f)/(1-rho)
        return(quan_tmp_u_r)
      }
    }
    
    optimize_interval = c(0,1)
    
    numCores = parallel::detectCores()
    res = mclapply(a,FUN=function(a.tmp){
      if(!is.null(grid) & !is.null(quan_rho_grid_inner)){
        grid_inner = as.numeric(names(quan_rho_grid_inner))
        ### refne interval for optimization by running the following over a grid (no 0 and 1)
        grid_inner_index = order((quan_rho_grid_inner - grid_inner*a.tmp)/(1-grid_inner))[1]
        optimize_interval = c(grid[grid_inner_index],grid[grid_inner_index+2])
      }
      min_quantile = optimize(quan_rho,interval=optimize_interval,df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective
      tail.prob = 0
      if(class(try(CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)))!="try-error") tail.prob = CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)$Qq
      if(tail.prob<=acc) tail.prob = pvalue.liu(Q=min_quantile,weight=lambda_r,method=ChisqApp)$pvalue
      (1-tail.prob) * dchisq(a.tmp,df=df_f)
    },
    mc.cores=min(numCores,max_core))
    res = as.numeric(unlist(res))
    return(res)
  }
  
  pvalue = 1 - integrate(integrand,lower=0,upper=qchisq(min_pvalue,df=df_f,lower.tail=FALSE),min_pvalue=min_pvalue,df_f=df_f,lambda_r=lambda_r,acc=acc,grid=quan_grid[[1]],quan_rho_grid_inner=quan_grid[[2]],max_core=max_core)$value
  if(pvalue<0) pvalue = abs(pvalue)
  return(pvalue)
}

pvalue_minimizePValue_liu = function(u_f,u_r,df_f,lambda_r,min_pvalue=NULL,acc=5e-10,lim=1e4,ChisqApp="4M"){
  if(is.null(min_pvalue)){
    p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
      pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                      lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
      if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,
                                          weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),
                                          method=ChisqApp)$pvalue
    }
    min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=df_f,lambda_r=lambda_r,u_f=u_f,u_r=u_r,acc=acc)$objective
    min_pvalue = min(c(min_pvalue_inside01,pvalue.r,pvalue.f))
  }
  
  integrand = function(a,min_pvalue,df_f,lambda_r,acc=5e-10){
    quan_rho = function(rho,df_f,lambda_r,p_rho_min,u_f){
      if(rho==1) stop("rho is not allowed to be 1.")
      else{
        if(ChisqApp=="4M") parameters = WeightSumChisq_mod(weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),df=rep(1,length(weight)))
        if(ChisqApp=="3M") parameters = WeightSumChisq(weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),df=rep(1,length(weight)))
        # mean.Q, SD.Q, df, noncentral, mean.X, SD.X
        quan_tmp = (qchisq(p_rho_min,df=parameters["df"],ncp=parameters["noncentral"],lower.tail=FALSE) - parameters["mean.X"]) / parameters["SD.X"] * parameters["SD.Q"] + parameters["mean.Q"] # Note: use upper tail for quantile calcualtion to avoid Inf when p_rho_min is extremely small.
        names(quan_tmp) = c()
        quan_tmp_u_r = (quan_tmp - rho*u_f)/(1-rho)
        return(quan_tmp_u_r)
      }
    }
    res = sapply(a,function(a.tmp){
      # min_quantile = optimize(quan_rho,interval=c(0,1-1e-3),df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective
      min_quantile_tmp = optimize(quan_rho,interval=c(0,1),df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp)$objective # try 2017-04-06
      min_quantile = min(quan_rho(rho=0,df_f=df_f,lambda_r=lambda_r,p_rho_min=min_pvalue,u_f=a.tmp),
                         min_quantile_tmp) # try 2017-04-06
      tail.prob = 0
      if(class(try(CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)))!="try-error") tail.prob = CompQuadForm::davies(q=min_quantile,lambda=lambda_r,acc=acc,lim=1/acc)$Qq
      if(tail.prob<=acc) tail.prob = pvalue.liu(Q=min_quantile,weight=lambda_r,method=ChisqApp)$pvalue
      (1-tail.prob) * dchisq(a.tmp,df=df_f)
    })
    return(res)
  }
  pvalue = 1 - integrate(integrand,lower=0,upper=qchisq(min_pvalue,df=df_f,lower.tail=FALSE),min_pvalue=min_pvalue,df_f=df_f,lambda_r=lambda_r,acc=acc)$value
  return(pvalue)
}

pvalue.liu = function(Q,weight,method="3M",df=rep(1,length(weight))){
  if(method=="4M") para = WeightSumChisq_mod(weight=weight,df=df)
  else para = WeightSumChisq(weight=weight,df=df)
  mean.Q = para[1]
  SD.Q = para[2]
  mean.X = para[5]
  SD.X = para[6]
  df.app = para[3]
  noncentral.app = para[4]
  Q.new = (Q-mean.Q)/SD.Q * SD.X + mean.X
  pvalue = pchisq(Q.new,df=df.app,ncp=noncentral.app,lower.tail=FALSE)
  return(list(pvalue=unname(pvalue),df=unname(df.app),ncp=unname(noncentral.app)))
}

Weight_ScrnStat = function(data,d=0,p,method){
  
  n = nrow(data)
  y = data[,1]
  e = data[,2]
  if(d!=0) {
    G = data[,(3+d):(3+d+p-1)]
    X = data[,3:(3+d-1)]
  }
  else G = data[,3:(3+p-1)]
  
  if(length(setdiff(c(0,1),e))==0)
  {
    cor.func = function(g,e){
      if(d!=0) margin.fit = summary(glm(e~g+X,family="binomial"))$coef
      else margin.fit = summary(glm(e~g,family="binomial"))$coef
      if("g" %in% rownames(margin.fit)) res = margin.fit[2,c(1,4,3)]
      else res = c(0,1,0)
      return(res)
    }
  }
  else
  {
    cor.func = function(g,e){
      if(d!=0) margin.fit = summary(lm(e~g+X))$coef
      else margin.fit = summary(lm(e~g))$coef
      if("g" %in% rownames(margin.fit)) res = margin.fit[2,c(1,4,3)]
      else res = c(0,1,0)
      return(res)
    }
  }
  cor.wt = apply(G,2,cor.func,e)
  
  marg.func = function(g,d){
    if(d!=0) summary(glm(y~g+X,family="binomial"))$coef[2,c(1,4,3)]
    else summary(glm(y~g,family="binomial"))$coef[2,c(1,4,3)]
  }
  marg.wt = apply(G,2,marg.func,d)
  
  if(method=="ScrnStat1"){
    wt = cor.wt
    wt[,cor.wt[2,]>=marg.wt[2,]] = marg.wt[,cor.wt[2,]>=marg.wt[2,]]
    wt.return = wt[3,]
  }
  else if(method=="ScrnStat2"){
    wt.sign = sign(cor.wt[3,])
    wt.return = wt.sign* (cor.wt[3,]^2+marg.wt[3,]^2)
  }
  else stop("method needs to be either ScrnStat1 or ScrnStat2")
  
  return(wt.return)
}

#R/SurvivalKMTests-v3-10102015.R
KMTest.surv = function(Z=NULL, U, Delta, X = NULL, gamma = NULL, kernel, npert = 1000){
  #=======================================================================================================================#
  #=======================================================================================================================#
  # Version 3: 10.10.2015                                                                                                 #
  # Changes in Version 3                                                                                                  #
  # Minor changes to increase efficiency, shd give identical results as V2                                                #
  #                                                                                                                       #
  # Version 2: 09.30.2013                                                                                                 #
  # Changes in Version 2:                                                                                                 #
  # doesn't check if Z, U, Delta, X etc are matrices etc                                                                  #
  # kernel must be an n x n kernel matrix that is supplied                                                                #
  # Z is not needed anymore, and removed the weights argument                                                             #
  #                                                                                                                       #
  # Notes: This function tests for association between a set of SNPS and a censored survival outcome                      #
  #        using the survival kernel machine.  The IBS and weighted kernels assume that the SNPs are in additive mode,    # 
  #        but can easily be adapted for SNPs in dominant mode by either defining new kernel functions or entering the    #
  #        nXn kernel matrix directly.                                                                                    #
  #        The function also assumes that survival times are right-censored.                                              # 
  #        The function also assumes that there are no missing data, individuals with missing data should be removed      #
  #        prior to analysis.                                                                                             #
  #                                                                                                                       #
  #                                                                                                                       #
  # :::INPUT:::                                                                                                           #
  #                                                                                                                       #
  # X:      is a nXR  matrix of relevant COVARIATES with each row as a different individual and each column               #
  #         as a separate covariate measurement. If no additional covariates are present, X can be left unspecified or    #
  #         left as NULL. Note that each column of X has to be a numerical variable, non-numerical variables have to be   #
  #         recoded appropriately before analysis.                                                                        #                         
  # Z:      is a nXS  matrix of the relevant SNPS with each row as a different individual and each column                 #
  #         as a separate SNP.  We asume the genotypes are in additive mode: 0,1,2 for non-carriers, heterozygotes        #
  #         and homozygous carriers, respectively.                                                                        #
  # U:      is a nX1 vector containing the observed times. Note: U=min(C,T) where C = censoring time, T = survival time   #
  # Delta:  is a nX1 vector containing the event indicator.                                                               #
  # gamma: Unless X = NULL, gamma has to be supplied.     gamma <- coxph(Surv(U,Delta)~X)$coef (?????)                    #                         
  # kernel: defines the kernel matrix.  This may be defined in several ways:                                              #
  #         kernel can be the nXn Kernel matrix OR                                                                        #
  #         kernel can be:                                                                                                #
  #                       "linear":  Linear Kernel                                                                        #
  #                       "IBS": IBS Kernel                                                                               #
  #                       "IBS.weighted": weighted IBS (weights must be specified a priori)                               #
  # weights: is a vector of length S*n of prespecified weights for the weighted IBS Kernel                                #
  #          See Kern.FUN.weightedIBS() to for information on how weights are defined                                     #
  # npert:   is the number of perturbations used to calculate p-value (default =1000)                                     #
  #                                                                                                                       #
  #                                                                                                                       #                                                                                                                                                                                                                                    
  # :::OUTPUT:::                                                                                                          #
  #                                                                                                                       # 
  # n : no. of samples                                                                                                    #
  # nsnps: no. of SNPs = S                                                                                                #   
  # p.value.ptb: single p-value for testing the null of no association based on resampling                                #
  # pvalue.chisq: single p-value for testing the null of no association based on Satterthwaite approximation              #
  #               Note that the Satterthwaite approximation uses resampling to estimate the scale parameter and df,       #
  #               and changes with each run. Use set.seed() to obtain identical p-values.                                 #
  # df.chisq: the degrees of freedom of the chi-square statistic                                                          #
  # score : the unscaled score statistic                                                                                  #
  # wptb: uncentered realizations of score statistic under null                                                           #                                                                                                                      
  # Citation:                                                                                                             #
  #                                                                                                                       #
  #                                                                                                                       #
  # If you use the survival kernel machine, please cite the following two papers:                                         #
  #                                                                                                                       #
  # Cai T, Tonini G and Lin X. 2011. Kernel machine approach to testing the significance of genomic pathways.             #
  # Biometrics 67: 975-86.                                                                                                #
  #                                                                                                                       #
  # Lin X, Cai T, Wu M, Zhou Q, Liu G, Christiani D, Lin X. 2011. Kernel Machine SNP-set Analysis for Censored            #
  # Survival Outcomes in Genome-wide Association Studies. Genetic Epidemiology 35: 620-31.                                #
  #                                                                                                                       #
  #                                                                                                                       #
  #                                                                                                                       #
  #=======================================================================================================================#
  #=======================================================================================================================#
	
      #if (!is.null(X)) {
   	  #   if (class(X)!= "matrix") stop("X is not a matrix")
	  #   if (nrow(X)!=nrow(Z)) stop("Dimensions of Z and X do not match")
      #}
  
      #if (class(Z)!= "matrix") stop("Z is not a matrix")
      #if (nrow(Z)!=length(U)) stop("Dimensions of Z and U do not match")
      #if (nrow(Z)!=length(Delta)) stop("Dimensions of Z and Delta do not match")


	# ====================================== #
	# Construct Kernel Matrix 
	# ====================================== #
	

		Kij = kernel


	nn = length(U); Gi = matrix(rnorm(npert*nn),nrow=nn)
	
	tl = U[Delta==1] # failure time
	stl = sort(tl); nl = length(tl) 
	if (is.null(X)) V = rep(1,nn) else V = as.vector(exp(X%*%gamma))
	pi0hat.stl = sum.I(stl,"<=", U, V) ## nl X 1
	dLam.stl <- 1/pi0hat.stl; Lam.stl <- cumsum(dLam.stl) ## nl X 1   
	# For each subject, tmpind is the number of the failure time is smaller than the time
    tmpind <- sum.I(U,">=",stl)+1; Lam.U <- c(0,Lam.stl)[tmpind] # sum{I(Xi>=stl)dLam.stl} n X 1
    Mhati <- Delta - Lam.U*V ## n X 1
    

    
    if (!is.null(X)){
    	q0 = ncol(X)
    	# S^(1), S^(1) S'^(1), S^(2)
    	pi1hat.stl <- sum.I(stl,"<=",U,X*V)
    	pi1hat2.stl <- matrix(0,length(stl),q0^2)
    	for (k in 1:length(stl)) {pi1hat2.stl[k,] <- as.vector(pi1hat.stl[k,]%*%t(pi1hat.stl[k,]))}
    	X2.mat <- matrix(0,nn,q0^2)
    	for (i in 1:nn){ X2.mat[i,] <- as.vector(X[i,]%*%t(X[i,]))}
    	pi2hat.stl = sum.I(stl,"<=",U,X2.mat*V)
    	Ahat = matrix(apply(pi2hat.stl/pi0hat.stl-pi1hat2.stl/pi0hat.stl^2,2,sum),q0,q0,byrow=T)/nn

	   #-----------------------------------------------
	   # change in V3 - save invAhat
       invAhat <- solve(Ahat)
	   #-----------------------------------------------
	   	    	
       W.gamma = W.gamma.t1 = W.gamma.t2 = W.gamma.t3 = matrix(0,nn,q0) 
       W.gamma.t1 = X*Mhati
       W.gamma.t2[Delta==1,] = pi1hat.stl[rank(tl),]/pi0hat.stl[rank(tl)]
       int.t3 = apply(pi1hat.stl/(pi0hat.stl)^2,2,cumsum)
       tmpind <- sum.I(U,">=",stl)+1
       W.gamma.t3 = rbind(rep(0,q0),int.t3)[tmpind,]*V
       #-----------------------------------------------
	   # change in V3 - save invAhat
	   W.gamma <- t( invAhat%*%t(W.gamma.t1-W.gamma.t2+W.gamma.t3) )
       #W.gamma = t(solve(Ahat)%*%t(W.gamma.t1-W.gamma.t2+W.gamma.t3))   
       #-----------------------------------------------
       X.tilde <- X*(Lam.U*V) 
     	}

	#-----------------------------------------------     	
    # change in V3 - moved this block below 	
    term1.ii <- c(0,cumsum(dLam.stl/pi0hat.stl))[tmpind]
    term2.ij <- pmin(matrix(U,nrow=nn,ncol=nn),matrix(U,nrow=nn,ncol=nn,byrow=T))
    term2.ij <- matrix(sum.I(c(term2.ij),">=",stl,dLam.stl/pi0hat.stl),ncol=nn)
    Mhati.Mhatl = Mhati%*%t(Mhati)
	#----------------------------------------------- 
	   
    score = 0; wptb = rep(0,npert)
    term1 <- sum(diag(Kij)*(diag(Mhati.Mhatl)- (Lam.U*V - term1.ii*V^2)))
    term2 <- sum((Kij-diag(diag(Kij)))*(term2.ij*(V%*%t(V))+Mhati.Mhatl))
    score <- term1+term2
        
    Keig<-eigen(Kij,symmetric=T); Keig.value<-pmax(Keig$values,0); Keig.vec <- Keig$vectors # n times m
    ri.mat1 = Keig.vec*Mhati
    omega.stl = sum.I(stl,"<=",U,Keig.vec*V)/pi0hat.stl #nt x m matrix 
    ri.mat2 = Keig.vec*0 # n x m matrix
    ri.mat2[Delta==1,] = omega.stl[rank(tl),] 
    ri.mat2 = ri.mat2 - sum.I(U, ">=", stl, omega.stl*dLam.stl)*V # n \times m     
    ri.mat = ri.mat1 - ri.mat2
    if (!is.null(X)){
    	ri.mat2.3 = t(pi1hat.stl%*%t(W.gamma)*dLam.stl)%*%omega.stl/nn
       ri.mat3 = W.gamma%*%t(X.tilde)%*%Keig.vec/nn
       ri.mat = ri.mat + ri.mat2.3 - ri.mat3
    }
    out = VTM(sqrt(Keig.value),nn)*ri.mat
    wptb = apply((t(out)%*%Gi)^2,2,sum)
 	
 	## p value obtained from perturbation
 	wptb.c = wptb - mean(wptb)
 	p.value.ptb = mean((wptb.c-score)>0)
 	
 	## p value obtained from chi-square approximation
 	wptb.mu = mean(wptb); wptb.var = var(wptb)
 	chisq.scale = wptb.var/wptb.mu/2
 	chisq.df = wptb.mu/chisq.scale
 	p.value.chisq = 1-pchisq((score+wptb.mu)/chisq.scale,df=chisq.df)
 	
 	return(list(n=nn,nsnps=ncol(Z), shat=score, sptb = wptb, pvalue.ptb=p.value.ptb,pvalue.Q=p.value.chisq, df.chisq=chisq.df))     
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
{
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    # for each distinct ordered failure time t[j], number of Xi < t[j]
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
    if (!is.null(Vi)) {
       ## if FUN contains '=', tmpind is the order of decending
        if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
        ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
        Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
        return(rbind(0,Vi)[pos+1,])
    } else return(pos)
}

Kern.FUN.Lin <- function(zz,weights=NULL)
{
    zz = as.matrix(zz);
    zz%*%t(zz)
}

# Kern.FUN.IBS modified April 28, 2009
Kern.FUN.IBS <- function(zz,weights=NULL)
  {
    # computes the average IBS between two individuals 
    # sum over IBS at all markers and divide by total no. of markers
    # only non-missing markers used (missing value coded as NA)
    # missing data is not included in both the numerator and denominator
    # max possible value is 1 if IBS=2 at ALL non-missing markers
    temp <- zz
    Ina <- 1*(is.na(temp)==F)
    zz[is.na(zz)] <- -9
    I0 = 1*(zz==0); I1 = 1*(zz==1); I2 = 1*(zz==2); 
    maxsnps <- Ina%*%t(Ina)
    (I0%*%t(I0)+I1%*%t(I1)+I2%*%t(I2)+0.5*I1%*%t(I0+I2)+0.5*(I0+I2)%*%t(I1))/max(maxsnps,0.00001)
  }

# Kern.FUN.weightedIBS modified June 28, 2010
Kern.FUN.weightedIBS <- function(zz, weights)
  {
    # computes the average IBS between two individuals 
    # weights is a vector of length p*n (not a matrix)
    # if there are S SNPs and n samples and we are weighting inversely by sqrt(MAF),
    # weights is of the form c(rep(MAF1, nsamples), rep(MAF2, nsamples), ... , rep(MAFp, nsamples))
    # sum over weighted IBS at all markers and divide by total no. of markers
    # input matrix CANNOT contain missing genotypes
    # max possible value is 1 if IBS=2 at ALL markers

    I0 = 1*(zz==0); I1 = 1*(zz==1); I2 = 1*(zz==2); 
    I0 = I0*(1/weights)^0.25; I1 = I1*(1/weights)^0.25; I2 = I2*(1/weights)^0.25;
    I0%*%t(I0)+I1%*%t(I1)+I2%*%t(I2)+0.5*I1%*%t(I0+I2)+0.5*(I0+I2)%*%t(I1)
  }


# Date: June 5, 2009
mafcall <- function(x){
	min(sum(x, na.rm=T), 2*sum(is.na(x)==F)-sum(x, na.rm=T))/(2*sum(is.na(x)==F)) 
}

# Date: June 5, 2009
checkpolymorphic <- function(x){
	length(unique(x))>1
}

VTM<-function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


# R/Kernel.R
# Sept 30, 2013
# Kernel.R in coxSKAT is modified from Kernel.R in SKAT (version 0.91).
# The only modification is to lskmTest.GetKernel()
# to add the option of linear and weighted linear kernels.
# The other functions in Kernel.R are essentially unchanged, except PACKAGE="coxSKAT" 
# is added whenever .C() is called.
# kernel_func.c is copied unchanged from (version 0.91).

K1_Help= function(x,y){
  # Helper function for 2 way interaction kernel
  p = length(x)
  a = x*y
  b = cumsum(a)
  return(sum(a[-1]*b[-p]))
}

call_Kernel_IBS<-function(Z,n,p){

	#Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_IBS",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)),PACKAGE="coxKM")[[4]]
	matrix(temp,nrow=n)
}

call_Kernel_IBS_Weight<-function(Z,n,p,weights){

	#Kernel_IBS_Weight(int * Z, int * pn, int * pp, int *UseGivenWeight ,  double * weight, double * Kernel)
	given_weight = 1;
	if( is.null(weights)){
		weights = rep(0,p);
		given_weight = 0;
	} else {
		# change!!
		weights<-weights^2;
	}
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_IBS_Weight",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.integer(given_weight),
	as.double(weights),as.double(as.vector(K)), PACKAGE="coxKM")[[6]]
	matrix(temp,nrow=n)
}

call_Kernel_2wayIX<-function(Z,n,p){

	#Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_2wayIX",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)),PACKAGE="coxKM")[[4]]
	matrix(temp,nrow=n)
}

lskmTest.GetKernel = function(Z, kernel, weights,n,m){

# Start of change -----------------------------------
    	if (kernel == "linear") {
      		K = Z%*%t(Z)    	
	}


    	if (kernel == "linear.weighted") {
                	Zweighted <- t(t(Z) * (weights))
			K = Zweighted %*% t(Zweighted)    
			rm(Zweighted)
	}

# End of change -----------------------------------

    	if (kernel == "quadratic") {
      		K = (Z%*%t(Z)+1)**2
    	}


	if (kernel == "IBS") {
      		K = call_Kernel_IBS(Z,n,m)
    	}
    	if (kernel == "IBS.weighted") {

      		K = call_Kernel_IBS_Weight(Z,n,m,weights)
    	}
  	if (kernel == "2wayIX") {
      		K = call_Kernel_2wayIX(Z,n,m)
    	}  
   	if (kernel == "IBS.weighted_OLD") {
      		#K = matrix(nrow = n, ncol = n)
      		if (is.null(weights)) {
        		qs = apply(Z, 2, mean)/(2)
        		weights = 1/sqrt(qs)
      		} else {
			weights<-weights^2
		}
      		K1 = matrix(nrow =n, ncol = n)
      		for (i in 1:n) {
        		K1[i,] = apply(abs(t(Z)-Z[i,])*weights,2, sum)
      		}
      		K= 1-(K1)/(2*sum(weights))
    	}

    	if (kernel == "IBS_OLD") {
      		K1=matrix(nrow=n,ncol=n)
      		for (i in 1:n) {
        		K1[i,] = apply(abs(t(Z)-Z[i,]),2, sum)
      		}
      		K = (2*m-K1)/(2*m)
    	}
   	if (kernel == "2wayIX_OLD") {
      		K = 1+Z%*%t(Z)
      		N1=  matrix(nrow = n, ncol = n)
      		for (i in 1:n){
        		for (j in i:n){
	    			N1[j,i] = N1[i,j] = K1_Help(Z[i,], Z[j,])
	  		}
      		}
      		K = K+N1
    	}
	return(K)

}


# R/Main.R
# Last Edited: Oct 2015

coxKM <- function(Z=NULL, U, Delta, X = NULL, gamma = NULL, kernel = "linear", weights = NULL, npert = 10^4, 
			npert.check=TRUE, npert.upper=10^8, npert.threshold=50,
			impute.method = "fixed", is_check_genotype=TRUE, is_dosage=FALSE, missing_cutoff=0.15, SetID=NULL){




	  #-------------------------------------
          # check dimensions when kernel matrix is not supplied
	  # if kernel matrix is supplied, Z is not required

	  if(class(kernel) != "matrix"){ 
	  	
	  	
	  #-------------------------------------
	  # only change in V0.3 is to move this under the "if(class(kernel) != "matrix"){ " statement
      # check that weights are supplied when kernel="IBS.weighted" or "linear.weighted"
	  if(is.null(weights)==TRUE & kernel=="IBS.weighted") stop("For weighted kernels, weights must be specified.")
      if(is.null(weights)==TRUE & kernel=="linear.weighted") stop("For weighted kernels, weights must be specified.")
	  #-------------------------------------
	  	
	  	if(class(Z)!= "matrix") stop("Z is not a matrix.")
        	if(nrow(Z)!=length(U)) stop("Dimensions of Z and U do not match")
        	if(nrow(Z)!=length(Delta)) stop("Dimensions of Z and Delta do not match")
	  

	 	if(!is.null(X)){
	   	   if(class(X)!= "matrix") stop("X is not a matrix.")
        	   if(nrow(X)!= nrow(Z)) stop("Dimensions of X and Z don't match.")
	  	}


	  	if(is.null(weights)==FALSE){
        		if(ncol(Z)!= length(weights)) stop("Dimensions of Z and weights don't match.")
			if(sum(weights<0)!=0) stop("weights have to be non-negative and cannot be missing.")
	  	}
	  }



	  #-------------------------------------
          # check dimensions when kernel matrix is supplied
	  # if kernel matrix is supplied, Z is not required

	  if(class(kernel) == "matrix"){ 
	  	if(nrow(kernel)!= ncol(kernel)) stop("kernel is not a square matrix.")
        	if(nrow(kernel)!=length(Delta)) stop("Dimensions of U and Delta do not match")
        	if(length(U)!=length(Delta)) stop("Dimensions of U and Delta do not match")
	  

	 	 if(!is.null(X)){
	   	   if(class(X)!= "matrix") stop("X is not a matrix.")
        		if(nrow(X)!= length(U)) stop("Dimensions of X and U don't match.")
	  	}
	  }



	  #-------------------------------------
	  # check for missing values in X, U, Delta

	  if(is.null(X)==FALSE){
        	if(sum(is.na(X))!= 0) stop("X cannot have any missing values.")
	  }

	  if(sum(is.na(U))!= 0) stop("U cannot have any missing values.")
	  if(sum(is.na(Delta))!= 0) stop("Delta cannot have any missing values.")



	  #-------------------------------------
	  # check that X doesn't have intercept
	  if(is.null(X)==FALSE){

		if(ncol(X)==1){
			if(checkpolymorphic(X)==FALSE) stop("X should not include intercept and must have more than one level.")	
		}else{
        		if(sum(apply(X, 2, checkpolymorphic))!= ncol(X)) stop("X should not include intercept and must have more than one level.")

		}
	  }



	#-------------------------------------
	# check Z and impute when kernel matrix is not supplied

	if(class(kernel) != "matrix"){ 
 	if(is_dosage ==TRUE){
		impute.method="fixed"
	}

	if(is_check_genotype==TRUE | is_dosage==TRUE){
		Z.out <- coxSKAT_MAIN_Check_Z(Z=Z, SetID=SetID, weights=weights, impute.method=impute.method,  missing_cutoff=missing_cutoff)
		Z <- as.matrix(Z.out$Z.test)
		weights <- as.vector(Z.out$weights.Z.test)
	}
	
	if(is.null(Z)==TRUE){ 

		if(is.null(SetID)){
			msg <- sprintf("The Z matrix has no SNPs." )
		} else {
			msg <- sprintf("In %s, the Z matrix has no SNPs.", SetID )
		}

		warning(msg,call.=FALSE)
	      return(list(p.value=NA, Q=NA, n.marker.test=NA, n.indiv=NA, df= NA))
	}



	}



	#-------------------------------------
	# Get the kernel matrix

	if (class(kernel) == "matrix") {
		kernel.matrix <- kernel
	}else{
		kernel.matrix <- lskmTest.GetKernel(Z, kernel, weights, n=nrow(Z), m=ncol(Z)) # n x n, has to be a matrix
	}



	#-------------------------------------
	# Run coxSKAT

	if(npert.check==FALSE){
		if(npert<=10^4){
			coxSKAT.out <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=npert)
			return(list(p.value=max(0.5/npert, coxSKAT.out$pvalue.ptb), Q=coxSKAT.out$shat, n.marker.test=ncol(Z), 
				n.indiv=coxSKAT.out$n, df= coxSKAT.out$df.chisq,p.chisq=coxSKAT.out$pvalue.Q))
		} else{
			n.iterations <- ceiling(npert/10^4)
			p.iterated <- c()
			for(ttt in n.iterations){
				ptemp <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=10^4)
				p.iterated <- c(p.iterated, ptemp$pvalue.ptb)
			}
			return(list(p.value=max(0.5/(10^4*ceiling(npert/10^4)), mean(p.iterated)), Q=ptemp$shat, n.marker.test=ncol(Z),
                        	n.indiv=ptemp$n, df= ptemp$df.chisq,p.chisq=coxSKAT.out$pvalue.Q))

	
		}	
	}else{

	       coxSKAT.out <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=npert)
		if(coxSKAT.out$pvalue.ptb<=npert.threshold/npert){

			n.iterations <- ceiling(npert.upper/10^4)
                        p.iterated <- c()
                        for(ttt in n.iterations){
                                ptemp <- KMTest.surv(Z=NULL, U=U, Delta=Delta, X=X, gamma=gamma, kernel=kernel.matrix, npert=10^4)
                                p.iterated <- c(p.iterated, ptemp$pvalue.ptb)
                        }

                 return(list(p.value=max(0.5/(n.iterations*10^4),mean(p.iterated)), Q=coxSKAT.out$shat, n.marker.test=ncol(Z),
                                n.indiv=coxSKAT.out$n, df= coxSKAT.out$df.chisq,p.chisq=coxSKAT.out$pvalue.Q))


		}else{
			 return(list(p.value=max(0.5/npert, coxSKAT.out$pvalue.ptb), Q=coxSKAT.out$shat, n.marker.test=ncol(Z),
                                n.indiv=coxSKAT.out$n, df= coxSKAT.out$df.chisq,p.chisq=coxSKAT.out$pvalue.Q))
		}
	}



}














#     modified from iSKAT_MAIN_Check_Z()
#	Check the Z, and do imputation
#     iSKAT_MAIN_Check_Z() modified significantly from SKAT_MAIN_Check_Z from V0.78
#
coxSKAT_MAIN_Check_Z <- function(Z, SetID, weights=NULL, impute.method,  missing_cutoff){

	# check.Z.error = 0 : no snps removed, but some snps possibly imputed
	# check.Z.error = 1 : all snps removed, returns NULL matrix for Z
	# check.Z.error = 2 : some snps removed, remainder snps may have been imputed

	check.Z.error <- 0
	n <- nrow(Z)
	##############################################
	# Recode Missing to NA

	IDX_MISS <- union(which(is.na(Z)), which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS] <- NA
	} 

	###################################################
	# Check missing rates and exclude any SNPs with missing rate > missing_cutoff
	# Also exclude non-polymorphic SNPs
	m <- ncol(Z)
	ID_INCLUDE_SNP <- NULL
	for(i in 1:m){
		missing.ratio <- length(which(is.na(Z[,i])))/n
		sd1 <- sd(Z[,i], na.rm=TRUE)
		if(missing.ratio < missing_cutoff && sd1 > 0){
			ID_INCLUDE_SNP <- c(ID_INCLUDE_SNP,i)
		}
	}
	
	if(length(ID_INCLUDE_SNP) == 0){

		if(is.null(SetID)){
			msg <- sprintf("ALL SNPs have either high missing rates or no-variation. ")
		} else {
			msg <- sprintf("In %s, ALL SNPs have either high missing rates or no-variation. ",SetID )
		}
		warning(msg, call.=FALSE)
		
		re <- list(Z.test=NULL, weights.Z.test=NULL, check.Z.error=1) 
		  

	} else{

 		if(m - length(ID_INCLUDE_SNP) > 0 ){

			if(is.null(SetID)){
				msg <- sprintf("%d SNPs with either high missing rates or no-variation are excluded!", m - length(ID_INCLUDE_SNP))
			} else {
				msg <- sprintf("In %s, %d SNPs with either high missing rates or no-variation are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
			}

			warning(msg, call.=FALSE)	
			Z <- as.matrix(Z[,ID_INCLUDE_SNP])
			if(is.null(weights)==FALSE) weights <- weights[ID_INCLUDE_SNP]
			check.Z.error <- 2
			IDX_MISS <- which(is.na(Z))
		}
	
	

		##########################################
		# Missing Imputation
		if(length(IDX_MISS) > 0){
			if(is.null(SetID)){
				msg <- sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
			} else {
				msg <- sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
			}

			warning(msg,call.=FALSE)
			Z <- Impute_coxSKAT(Z,impute.method)
		} 
		re <- list(Z.test=Z, weights.Z.test=weights, check.Z.error=check.Z.error)
	
	}

	return(re)
}




# copied without modifcation from SKAT V0.78, renamed Impute_coxSKAT()
# Simple Imputation
# Z : an n x p genotype matrix with n samples and p SNPs
# Missing has to be NA: a missing genotype value.

Impute_coxSKAT <-function(Z, impute.method){
	
	p <- dim(Z)[2]

	if(impute.method =="random"){
		for(i in 1:p){
			IDX <- which(is.na(Z[,i]))
			if(length(IDX) > 0){
				maf1 <- mean(Z[-IDX,i])/2
				Z[IDX,i] <- rbinom(length(IDX),2,maf1)
			}
		}
	} else if(impute.method =="fixed"){
		for(i in 1:p){
			IDX<-which(is.na(Z[,i]))
			if(length(IDX) > 0){
				maf1 <- mean(Z[-IDX,i])/2
				Z[IDX,i] <- 2 * maf1
			}
		}
	} else {
		stop("Error: Imputation method shoud be either \"fixed\" or \"random\"! ")
	}

	return(Z)
}

  if(ncol(M)==0) stop("There must be at least 1 varant in the set to perform IECH.")
  if(ncol(M)==1) stop("Variance component test requires ncol(M)>=2.")
  
  
  n = nrow(data)
  EigenvalueD.r=rep(1,times=ncol(M))


 pvalue = c(pvalue.Burden=NA,pvalue.KM=NA,pvalue.IEHC_optim=NA,pvalue.IEHC_adapt=NA,pvalue.IEHC_Fisher=NA,pvalue.ACAT=NA)

  
  y_surv = Surv(data$time,data$status)
  fit0 = summary(coxph(y_surv ~ X + G))
  pvalue.Burden = fit0$coefficients[dim(X)[2]+1,5]

  Gamma0 <- coxph(Surv(data$time, data$status)~X)$coef
  pvalue.KM<-coxKM(Z=M, U=data$time, Delta=data$status, X,gamma=Gamma0,kernel="linear",is_check_genotype=FALSE)$`p.chisq`                     
  Gamma1 <- coxph(Surv(data$time, data$status)~cbind(X,G))$coef
  pvalue.chisq<-coxKM(Z=M, U=data$time, Delta=data$status, X=cbind(X,G),gamma=Gamma1,kernel="linear",is_check_genotype=FALSE)$`p.chisq`
  df<-coxKM(Z=M, U=data$time, Delta=data$status, X=cbind(X,G),gamma=Gamma1,kernel="linear",is_check_genotype=FALSE)$`df`
  Q.f=qchisq(pvalue.Burden, 1,  lower.tail = FALSE)
  Q.r=qchisq(pvalue.chisq, df,  lower.tail = FALSE)

      if(combinations_return & (("All" %in% combination_preference) | ("optim" %in% combination_preference))){
        # optim method - Calculate the optimal weight for linear combinations of Q.f and Q.r by minimizing the pvalues
        p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
          pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                          lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
          # Note on 2017/2/6: previous version doesn't include acc and lim input here and it causes the returned pvalue >1 in some simulation runs. Caused error in the following calculation of quantiles based on returned pvalues. acc and lim are added to control this.
          if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),method=chisq_app)$pvalue
          return(pvalue)
        }
        rho_min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=1,lambda_r=EigenvalueD.r,u_f=Q.f,u_r=Q.r,acc=acc)$minimum
        min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=1,lambda_r=EigenvalueD.r,u_f=Q.f,u_r=Q.r,acc=acc)$objective
        ### Note: optimize function does NOT include the comparison at the endpoints of the interval. Add comparison with p_values at 0 and 1 manually below
        if(min_pvalue_inside01>pvalue.chisq | min_pvalue_inside01>pvalue.Burden){
          rho_min_pvalue = ifelse(pvalue.chisq<pvalue.Burden,0,1)
          min_pvalue = ifelse(pvalue.chisq<pvalue.Burden,pvalue.chisq,pvalue.Burden)
        }
        else{
          rho_min_pvalue = rho_min_pvalue_inside01
          min_pvalue = min_pvalue_inside01
        }

        Q.Optim = Q.f*rho_min_pvalue + Q.r*(1-rho_min_pvalue)

  
          acc_supp = max(min(min_pvalue,0.01)*0.1,acc)
          lim_supp = as.integer(1/acc_supp)
 

        if(-log10(min_pvalue)>accurate_app_threshold){
          quan_grid = quan_rho_calculator(grid=seq(0,1,by=0.05),df_f=1,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc_supp,lim=lim_supp,max_core=max_core)
          pvalue.IEHC_optim = pvalue_minimizePValue_davies(u_f=Q.f,u_r=Q.r,df_f=1,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc,ChisqApp=chisq_app,max_core=max_core)
        }
        else pvalue.IEHC_optim = pvalue_minimizePValue_liu(u_f=Q.f,u_r=Q.r,df_f=1,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc,ChisqApp=chisq_app)
        pvalue["pvalue.IEHC_optim"] = pvalue.IEHC_optim
        #stat["stat.oMiST"] = Q.Optim
        #rho["rho.oMiST.f"] = rho_min_pvalue
        #rho["rho.oMiST.r"] = 1-rho_min_pvalue
      }

      if(combinations_return & (("All" %in% combination_preference) | ("adapt" %in% combination_preference))){
        # Adaptive weighted combination
        Z.f = -2*log(pvalue.Burden)
        if(is.null(log(pvalue.chisq))) {Z.r = -2*log(1e-16)}
        else {Z.r = -2*log(pvalue.chisq)}
        esti.weight.f = Z.f/sqrt(Z.f^2+Z.r^2)
        esti.weight.r = Z.r/sqrt(Z.f^2+Z.r^2)
        Q.Adapt = sqrt(Z.f^2+Z.r^2) # [1,] is added to avoid operation problems. It guarantees it's a value.
        integrand = function(y){
          pchisq(sqrt(Q.Adapt^2-y^2),df=2,lower.tail=FALSE)*dchisq(y,df=2)
        }
        pvalue.IEHC_adapt = pchisq(Q.Adapt,df=2,lower.tail=FALSE)+integrate(integrand,lower=0,upper=Q.Adapt)$value
        pvalue["pvalue.IEHC_adapt"] = pvalue.IEHC_adapt
        #stat["stat.aMiST"] = Q.Adapt
        #rho["rho.aMiST.f"] = esti.weight.f
        #rho["rho.aMiST.r"] = esti.weight.r
      }

      if(combinations_return & (("All" %in% combination_preference) | ("Fisher" %in% combination_preference))){
        # Fisher's combination
        Q.Fisher = -2*log(pvalue.Burden)-2*log(pvalue.chisq)
        pvalue.IEHC_Fisher = pchisq(as.numeric(Q.Fisher),df=4,lower.tail=FALSE)
        pvalue["pvalue.IEHC_Fisher"] = pvalue.IEHC_Fisher
        #stat["stat.fMiST"] = Q.Fisher
      }
      if(combinations_return & (("All" %in% combination_preference) | ("Burden" %in% combination_preference))){
        pvalue["pvalue.Burden"] = pvalue.Burden
      }
       if(combinations_return & (("All" %in% combination_preference) | ("KM" %in% combination_preference))){
        pvalue["pvalue.KM"] = pvalue.KM
      }   
       if(combinations_return & (("All" %in% combination_preference) | ("ACAT" %in% combination_preference))){
        Q.Fisher = -2*log(pvalue.Burden)-2*log(pvalue.chisq)
        pvalue.IEHC_Fisher = pchisq(as.numeric(Q.Fisher),df=4,lower.tail=FALSE)
	
	 Z.f = -2*log(pvalue.Burden)
        if(is.null(log(pvalue.chisq))) {Z.r = -2*log(1e-16)}
        else {Z.r = -2*log(pvalue.chisq)}
        esti.weight.f = Z.f/sqrt(Z.f^2+Z.r^2)
        esti.weight.r = Z.r/sqrt(Z.f^2+Z.r^2)
        Q.Adapt = sqrt(Z.f^2+Z.r^2) # [1,] is added to avoid operation problems. It guarantees it's a value.
        integrand = function(y){
          pchisq(sqrt(Q.Adapt^2-y^2),df=2,lower.tail=FALSE)*dchisq(y,df=2)
        }
        pvalue.IEHC_adapt = pchisq(Q.Adapt,df=2,lower.tail=FALSE)+integrate(integrand,lower=0,upper=Q.Adapt)$value
        
	
p_rho = function(rho,df_f,lambda_r,u_f,u_r,acc=5e-10){
          pvalue = CompQuadForm::davies(q=rho*u_f+(1-rho)*u_r,
                          lambda=c(rho*rep(1,df_f),(1-rho)*lambda_r),acc=acc,lim=1/acc)$Qq
          # Note on 2017/2/6: previous version doesn't include acc and lim input here and it causes the returned pvalue >1 in some simulation runs. Caused error in the following calculation of quantiles based on returned pvalues. acc and lim are added to control this.
          if(pvalue<=acc) pvalue = pvalue.liu(Q=rho*u_f+(1-rho)*u_r,weight=c(rho*rep(1,df_f),(1-rho)*lambda_r),method=chisq_app)$pvalue
          return(pvalue)
        }
        rho_min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=1,lambda_r=EigenvalueD.r,u_f=Q.f,u_r=Q.r,acc=acc)$minimum
        min_pvalue_inside01 = optimize(p_rho,interval=c(0,1),maximum=FALSE,df_f=1,lambda_r=EigenvalueD.r,u_f=Q.f,u_r=Q.r,acc=acc)$objective
        ### Note: optimize function does NOT include the comparison at the endpoints of the interval. Add comparison with p_values at 0 and 1 manually below
        if(min_pvalue_inside01>pvalue.chisq | min_pvalue_inside01>pvalue.Burden){
          rho_min_pvalue = ifelse(pvalue.chisq<pvalue.Burden,0,1)
          min_pvalue = ifelse(pvalue.chisq<pvalue.Burden,pvalue.chisq,pvalue.Burden)
        }
        else{
          rho_min_pvalue = rho_min_pvalue_inside01
          min_pvalue = min_pvalue_inside01
        }

        Q.Optim = Q.f*rho_min_pvalue + Q.r*(1-rho_min_pvalue)

  
          acc_supp = max(min(min_pvalue,0.01)*0.1,acc)
          lim_supp = as.integer(1/acc_supp)
 

        if(-log10(min_pvalue)>accurate_app_threshold){
          quan_grid = quan_rho_calculator(grid=seq(0,1,by=0.05),df_f=1,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc_supp,lim=lim_supp,max_core=max_core)
          pvalue.IEHC_optim = pvalue_minimizePValue_davies(u_f=Q.f,u_r=Q.r,df_f=1,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc,ChisqApp=chisq_app,max_core=max_core)
        }
        else pvalue.IEHC_optim = pvalue_minimizePValue_liu(u_f=Q.f,u_r=Q.r,df_f=1,lambda_r=EigenvalueD.r,min_pvalue=min_pvalue,acc=acc,ChisqApp=chisq_app)
        

	px=c(pvalue.Burden,pvalue.KM,pvalue.IEHC_Fisher,pvalue.IEHC_adapt,pvalue.IEHC_optim)
	pxx = na.omit(as.numeric(paste(px)))
	index = which(pxx==1)
	if (length(index) > 0) pxx[index] = 0.9999
	index = which(pxx<1e-30)
	if (length(index) > 0) pxx[index] = 1e-30
	pvalue.ACAT = ACAT(pxx)
	

	pvalue["pvalue.ACAT"] = pvalue.ACAT
      }   

 

  return(list(pvalue=pvalue))
}

