adjust_quantile <-
  function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues,method=0){
    
    ## input_pvalues is the 2-column matrix storing the two sets of p-values
    ## alpha00, alpha01, alpha10 are the proportions of three nulls, 
    ## method=0 corresponding to the approximation used in section 2.2-2.3 of the paper;
    ## method=1 corresponding to the exact method used in section 2.4 of the paper
    
    
    #library(fdrtool)
    nmed <- nrow(input_pvalues)  
    
    ## first compute the quantiles using the approximation method
    
    pnull <- rep(0,nmed) 
    for (i in 1:nmed) {
      c <- (-i/nmed)
      b <- alpha10+alpha01
      a <- alpha00
      pnull[i] <- (-b+sqrt(b^2-4*a*c))/(2*a)
    }  
    
    if (method==1) {
      cdf12 <- input_pvalues
      orderp1 <- input_pvalues[order(input_pvalues[,1]),1]
      orderp2 <- input_pvalues[order(input_pvalues[,2]),2]
      
      xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
      yy1 <- c(0,seq(1,nmed,by=1)/nmed)
      gfit1<- gcmlcm(xx1,yy1,type="lcm")
      xknots1 <- gfit1$x.knots[-1]
      Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)
      
      xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
      yy2 <- c(0,seq(1,nmed,by=1)/nmed)
      gfit2<- gcmlcm(xx2,yy2,type="lcm")
      xknots2 <- gfit2$x.knots[-1]
      Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)
      
      if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
      if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))
      
      gcdf1 <- orderp1
      gcdf2 <- orderp2
      
      orderq1 <- orderp1
      orderq2 <- orderp2
      
      difff <- 1
      ite <- 1
      
      while(abs(difff)>1e-6 & ite<10) {
        #cat(ite,"..")
        for (i in 1:length(xknots1)) {
          if (i==1) {
            gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]] 
          } else {   
            if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
              temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] 
              gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
            }
          }
        }
        
        for (i in 1:length(xknots2)) {
          if (i==1) {
            gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]] 
          } else {   
            if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
              temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] 
              gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
            } 
          }
        }
        
        cdf12[,1] <- ifelse(gcdf1>1,1,gcdf1)
        cdf12[,2] <- ifelse(gcdf2>1,1,gcdf2)
        
        pnull <- rep(0,nmed) 
        for (i in 1:nmed) {
          c <- (-i/nmed)
          b <- alpha10*cdf12[i,1]+alpha01*cdf12[i,2]
          a <- alpha00
          pnull[i] <- (-b+sqrt(b^2-4*a*c))/(2*a)
        }  
        
        difff <- max(max(orderq1-pnull),max(orderq2-pnull))
        
        orderq1 <- pnull
        orderq2 <- pnull
        ite <- ite+1
      }
    }
    return(pnull)
  }
