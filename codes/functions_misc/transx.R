transx <- function(x,tcode,lag){
  if(!is.matrix(x)) x<-as.matrix(x)
  # transform x series
  # -- tcodes:
  #       1 Level
  #       2 First Difference
  #       3 Second Difference
  #       4 Log-Level
  #       5 Log-First-Difference
  #       6 Log-Second-Difference
  #       7 Detrend Log Using 1-sided HP detrending for Monthly data
  #       8 Detrend Log Using 1-sided HP detrending for Quarterly data
  #       9 (1-L)/(1-L^12)
  #      10 Annualizing Monthly Data
  #      11 Annualizing Quarterly Data
  #
  #       Translated from the Gauss procs of Stock&Watson(2005),'Implications of
  #       dynamic factor models for VAR analysis'
  #       Dimitris Korobilis, June 2007
  #       
  #       Translated from Matlab Code by Dimitris Korobilis to R
  #       Maximilian Böck, November 2019
  
  small   = 1e-06
  relvarm = .00000075
  relvarq = .000625
  
  n <- nrow(x)
  y <- matrix(NA,n,1)
  
  if(tcode==1){
    y<-x
  }else if(tcode==2){
    y[(lag+1):n,] <- diff(x, lag=lag, differences=1)
  }else if(tcode==3){
    y[(lag*2+1):n,] <- x[(lag*2+1):n,] - 2*x[(lag+1):(n-lag),] + x[1:(n-lag*2),]
  }else if(tcode==4||tcode==5||tcode==6||tcode==7||tcode==8||tcode==9||tcode==10||tcode==11||tcode==12){
    if(min(x) < small){
      y[1:n,]<-NA
    }else{
      if(tcode==4) y[1:n,]         <- log(x)*100
      if(tcode==5) y[(lag+1):n,]   <- diff(log(x), lag=lag, differences=1)*100
      if(tcode==6) y[(2*lag+1):n,] <- diff(log(x), lag=lag, differences=2)*100
      if(tcode==7) y[1:n,] <- detrend(log(x),relvarm)[["xc"]]
      if(tcode==8) y[1:n,] <- detrend(log(x),relvarq)[["xc"]]
      if(tcode==9) y[14:n,] <- log(x)[14:n,]-log(x)[13:(n-1),]-log(x)[2:(n-12)]+log(x)[1:(n-13),]
      if(tcode==10) y[(lag+1):n,] = ((x[(lag+1):n,]/x[1:(n-lag),])^12-1)*100
      if(tcode==11) y[(lag+1):n,] = ((x[(lag+1):n,]/x[1:(n-lag),])^4-1)*100
      if(tcode==12) y[1:n,] <- log(100*x)
    }
  }
  return(y)
}
