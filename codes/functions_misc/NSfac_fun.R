NSfac_fun <- function(yld.data, maturities, lambda1 = 0.0609){
  #three Nelson-Siegel factors
  NS.lambda <- as.data.frame(matrix(NA,length(maturities), 3)) 
  dimnames(NS.lambda) <- list(paste0("Y_", maturities/12), paste0("ILS_", c("LEV", "SLP", "CURV")))
  NS.fac <- as.data.frame(matrix(NA,nrow(yld.data), 3))
  dimnames(NS.fac) <- list(as.character(time(yld.data)), paste0("ILS_", c("LEV", "SLP", "CURV")))
  
  for (jj in 1:length(maturities)){
    s11 <- (1-exp(-maturities[[jj]]*lambda1))/(maturities[[jj]] * lambda1)
    c11 <- s11-exp(-maturities[[jj]]*lambda1)
    NS.lambda[jj,] <- c(1, s11, c11)
  }
  # Obtain OLS estimates for Nelson-Siegel factors
  for(tt in 1:nrow(yld.data)) NS.fac[tt,] <- t(solve(crossprod(as.matrix(NS.lambda))))%*%crossprod(as.matrix(NS.lambda),as.matrix(as.numeric(yld.data[tt,])))
  
  return(list(NS.fac=NS.fac, NS.lambda=NS.lambda))
}