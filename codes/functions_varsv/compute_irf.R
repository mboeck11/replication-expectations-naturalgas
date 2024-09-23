compute_irf <- function(run, nhor, scaling=NULL, ident="chol"){
  
  # dimensions
  n     = ncol(run$Y)
  k     = ncol(run$X)
  bigT  = nrow(run$Y)
  plag  = run$args$plag
  SV    = run$args$SV
  cons  = run$args$cons
  trend = run$args$trend
  quadratictrend = run$args$quadratictrend
  dets = ifelse(cons,1,0) + ifelse(trend,1,0) + ifelse(quadratictrend,1,0)
  
  # other
  thindraws = run$args$thindraws
  varNames  = colnames(run$Y)
  
  irf_store  = array(NA_real_, c(thindraws, n, n, nhor),
                     dimnames=list(NULL, varNames, varNames, seq(0,nhor-1)))
  for(irep in 1:thindraws){
    
    if(irep%%50==0) cat(paste0("IRF Computation: ",irep," / ", thindraws,".\n"))
    
    temp    <- gen_compMat(A=run$A[irep,1:(n*plag),], n=n, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    
    if(max(abs(Re(eigen(compMat)$values))) > 1)
      next
    
    shock <- t(chol(run$SIGMA[irep,,]))
    if(!is.null(scaling)) shock <- shock%*%diag(1/diag(shock))*scaling
    
    impresp1 <- array(NA_real_, c(n, n, nhor))
    impresp1[,1:n,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp1[,1:n,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    
    # save stuff
    irf_store[irep,,,]  = impresp1
  }
  idx <- which(!is.na(irf_store[,1,1,1]))
  thindraws  = length(idx)
  irf_store  = irf_store[idx,,,,drop=FALSE]
  
  return(irf_store=irf_store)
}
