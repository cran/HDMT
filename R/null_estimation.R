null_estimation <- function(input_pvalues,lambda=0.5) {
  
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
 
  if (ks.test(input_pvalues[,1],"punif",0,1,alternative="greater")$p>0.05) alpha1 <- 1 else   alpha1 <- min(frac1[pcut==lambda],1)  
  if (ks.test(input_pvalues[,2],"punif",0,1,alternative="greater")$p>0.05) alpha2 <- 1 else   alpha2 <- min(frac2[pcut==lambda],1)
  

  if (alpha00==1) {
    alpha01 <- 0
    alpha10 <- 0
    alpha11 <- 0
  } else {    
    if (alpha1==1  & alpha2==1) {
      alpha01 <- 0
      alpha10 <- 0
      alpha11 <- 0
      alpha00 <- 1
    }  
    
    if (alpha1==1  & alpha2!=1) {
      alpha10 <- 0
      alpha11 <- 0
      alpha01 <- alpha1-alpha00
      alpha01 <- max(0,alpha01)
      alpha00 <- 1-alpha01
    }  
    
    if (alpha1!=1  & alpha2==1) {
      alpha01 <- 0
      alpha11 <- 0
      alpha10 <- alpha2-alpha00
      alpha10 <- max(0,alpha10)
      alpha00 <- 1-alpha10
    }  
    
    if (alpha1!=1  & alpha2!=1) {
      alpha10 <- alpha2-alpha00
      alpha10 <- max(0,alpha10)
      alpha01 <- alpha1-alpha00
      alpha01 <- max(0,alpha01)
      
      if ((1-alpha00-alpha01-alpha10)<0) {
        alpha11 <- 0
        alpha10 <- 1- alpha1
        alpha01 <- 1- alpha2
        alpha00 <- 1- alpha10 - alpha01
      }  else {
        alpha11 <-  1-alpha00-alpha01-alpha10
      }  
    }  
  }
  alpha.null <- list(alpha10=alpha10,alpha01=alpha01,alpha00=alpha00,alpha1=alpha1,alpha2=alpha2)
  return(alpha.null)
}
