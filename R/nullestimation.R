nullestimation <-function(input_pvalues,lambda=0.5) {
  
  ## input_pvalues is a matrix with 2 columns of p-values, the first column is p-value for exposure-mediator association, the second column is p-value for mediator-outcome association adjusted for exposure
  ## lambda is the threshold for pi_{00} estimation, default 0.5
  #check input
  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")
  
  pcut <- seq(0.1,0.8,0.1) 
  frac1 <- rep(0,8)
  frac2 <- rep(0,8)
  frac12<- rep(0,8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[,1]>=pcut[i])/(1-pcut[i])
    frac2[i] <- mean(input_pvalues[,2]>=pcut[i])/(1-pcut[i]) 
    frac12[i]<- mean(input_pvalues[,2]>=pcut[i] & input_pvalues[,1]>=pcut[i])/(1-pcut[i])^2
  }  
  
  ## use the median estimates for pi00 ##
  
  alpha00 <- min(frac12[pcut==lambda],1)
  
  ## alpha1 is the proportion of nulls for first p-value 
  ## alpha2 is the proportion of nulls for second p-value 
  #library(cp4p)
  if (ks.test(input_pvalues[,1],"punif",0,1)$p>0.05) alpha1 <- 1 else     alpha1 <- estim.pi0(p=input_pvalues[,1], pi0.method="slim")$pi0.SLIM
  if (ks.test(input_pvalues[,2],"punif",0,1)$p>0.05) alpha2 <- 1 else     alpha2 <- estim.pi0(p=input_pvalues[,2], pi0.method="slim")$pi0.SLIM
  
  if (alpha1 == 1 | alpha2==1) {
    alpha11 <- 0
    alpha10 <- 1-alpha1
    alpha01 <- 1-alpha2
    if (alpha1==1 & alpha2==1) alpha00 <- 1 else {
      if (sum(alpha00+alpha10+alpha01)>1) {
        err <- sum(alpha00+alpha10+alpha01)-1
        alpha00 <- alpha00- err 
      }  else {
        err <- 1- sum(alpha00+alpha10+alpha01)
        if (alpha10 !=0) alpha10 <- alpha10 + err
        if (alpha01 !=0) alpha01 <- alpha01 + err
      }      
    }
  } else {
    
    alpha11 <- max(1-alpha1+1-alpha2+alpha00-1,0)
    if (alpha11 ==0) {
      alpha10 <- 1-alpha1 
      alpha01 <- 1-alpha2 
      err <- 1-(alpha00+alpha10+alpha01)
      alpha00 <- alpha00 + err
      
    }  else {
      if (alpha11> (1-alpha1)| alpha11> (1-alpha2)) {
        alpha11 <- 0
        alpha10 <- 1-alpha1
        alpha01 <- 1-alpha2
        
        err <- 1-(alpha00+alpha10+alpha01)
        alpha00 <- alpha00 + err
        
      } else {
        alpha10 <- max(1-alpha1-alpha11,0)
        alpha01 <- max(1-alpha2-alpha11,0)
        alpha00 <- 1-alpha10-alpha01-alpha11
      }
    }
  }
  alpha.null <- list(alpha10=alpha10,alpha01=alpha01,alpha00=alpha00,alpha1=alpha1,alpha2=alpha2)
  return(alpha.null)
}

