compute_hd <- function(run, compute_mean=TRUE){
  # get data 
  Y = run$Y
  X = run$X
  
  # get coefs
  A     = run$A
  Sigma = run$SIGMA
  
  thindraws = run$args$thindraws
  
  # dimensions
  n      = ncol(Y)
  k      = ncol(X)
  bigT   = nrow(Y)
  plag   = run$args$plag
  K      = n*plag
  SV     = run$args$SV
  cons   = run$args$cons
  trend  = run$args$trend
  qtrend = run$args$qtrend
  exo    = ifelse(is.null(run$args$Ex),FALSE,TRUE)
  Kex    = ifelse(!is.null(run$args$Ex),ncol(run$args$Ex),0)
  dets = ifelse(cons,1,0) + ifelse(trend,1,0) + ifelse(qtrend,1,0) + Kex
  
  # identification (Cholesky)
  R    = diag(n) # for sign-restrictions when we need it
  
  if(compute_mean){
    
    # kick non-stationary draws
    idx = NULL
    for(irep in 1:thindraws){
      temp = gen_compMat(A=run$A[irep,1:(n*plag),], n=n, p=plag)
      Cm   = temp$Cm
      Jm   = temp$Jm
      if(max(abs(Re(eigen(Cm)$values))) > 1)
        next
      
      idx = c(idx,irep)
    }
    
    Amean     = apply(run$A[idx,,], c(2,3), mean)
    Sigmamean = apply(run$SIGMA[idx,,], c(2,3), mean)
    
    temp = gen_compMat(A=Amean[1:(n*plag),], n=n, p=plag)
    Cm   = temp$Cm
    Jm   = temp$Jm
    
    Sigchol    = t(chol(Sigmamean))
    SigcholR   = Sigchol %*% R
    strMat     = solve(SigcholR)
    
    # structural errors
    eps = strMat %*% t(Y - X%*%run$A[irep,,])
    
    # contribution of each shock
    invA_big = matrix(0, K, n)
    invA_big[1:n,] = strMat
    HDshock_big = array(0, c(K, bigT+1, n))
    HDshock     = array(0, c(n, bigT+1, n))
    for(nn in 1:n){
      eps_big = matrix(0, n, bigT+1)
      eps_big[nn,2:(bigT+1)] = eps[nn,]
      for(tt in 2:(bigT+1)){
        HDshock_big[,tt,nn] = invA_big %*% eps_big[,tt] + Cm %*% HDshock_big[,(tt-1),nn]
        HDshock[,tt,nn]     = t(Jm) %*% HDshock_big[,tt,nn]
      }
    }
    hd_shock = aperm(HDshock[,2:(bigT+1),], c(2,3,1))
    
    # contribution of initial value
    HDinit_big = matrix(0, n*plag, bigT+1)
    HDinit     = matrix(0, n, bigT+1)
    HDinit_big[,1] = t(X[1,1:n*plag,drop=FALSE])
    HDinit[,1]     = t(Jm) %*% HDinit_big[,1]
    for(tt in 2:(bigT+1)){
      HDinit_big[,tt] = Cm %*% HDinit_big[,tt-1]
      HDinit[,tt]     = t(Jm) %*% HDinit_big[,tt]
    }
    hd_init  = t(HDinit[,2:(bigT+1)])
    
    # contribution of constant
    HDcons = matrix(0, n, bigT+1)
    if(cons){
      HDcons_big = matrix(0, n*plag, bigT+1)
      CC = matrix(0, K, 1)
      CC[1:n,] = A[irep,"cons",]
      for(tt in 2:(bigT+1)){
        HDcons_big[,tt] = CC + Cm %*% HDcons_big[,tt-1]
        HDcons[,tt]     = t(Jm) %*% HDcons_big[,tt]
      }
      hd_cons = t(HDcons[,2:(bigT+1)])
    }else{
      hd_cons = NULL
    }
    
    # contribution of trend
    HDtrend = matrix(0, n, bigT+1)
    if(trend){
      HDtrend_big = matrix(0, K, bigT+1)
      TT = matrix(0, K, 1)
      TT[1:n,] = A[irep,"trend",]
      for(tt in 2:(bigT+1)){
        HDtrend_big[,tt] = TT*(tt-1) + Cm %*% HDtrend_big[,tt-1]
        HDtrend[,tt]     = t(Jm) %*% HDtrend_big[,tt]
      }
      hd_trend = t(HDtrend[,2:(bigT+1)])
    }else{
      hd_trend = NULL
    }
    
    # contribution of quadratic trend
    HDqtrend = matrix(0, n, bigT+1)
    if(qtrend){
      HDqtrend_big = matrix(0, K, bigT+1)
      TT2 = matrix(0, K, 1)
      TT2[1:n,] = A[irep,"qtrend",]
      for(tt in 2:(bigT+1)){
        HDqtrend_big[,tt] = TT2*(tt-1)^2 + Cm %*% HDqtrend_big[,tt-1]
        HDqtrend[,tt]     = t(Jm) %*% HDqtrend_big[,tt]
      }
      hd_qtrend = t(HDqtrend[,2:(bigT+1)])
    }else{
      hd_qtrend = NULL
    }
    
    # contributio  of exogenous regressors
    HDexo = matrix(0, n, bigT+1)
    if(exo){
      HDexo_big = matrix(0, K, bigT+1)
      EXO = matrix(0, K, Kex)
      EXO[1:n,] = A[irep,grepl("Ex",dimnames(A)),]
      for(tt in 2:(bigT+1)){
        HDexo_big[,tt] = t(EXO %*% X[,grepl("Ex",colnames(X))]) + Cm %*% HDexo_big[,tt-1]
        HDexo[,tt]     = t(Jm) %*% HDexo_big[,tt]
      }
      hd_exo = t(HDexo[,2:(bigT+1)])
    }else{
      hd_exo = NULL
    }
    
    return(list(hd_shock=hd_shock, hd_init=hd_init, hd_cons=hd_cons,
                hd_trend=hd_trend, hd_qtrend=hd_qtrend, hd_exo=hd_exo))
  }else{
    hd_store     = array(NA_real_, c(thindraws, bigT, n, n))
    init_store   = array(NA_real_, c(thindraws, bigT, n))
    if(cons)   cons_store   = array(NA_real_, c(thindraws, bigT, n)) else cons_store = NULL
    if(trend)  trend_store  = array(NA_real_, c(thindraws, bigT, n)) else trend_store = NULL
    if(qtrend) qtrend_store = array(NA_real_, c(thindraws, bigT, n)) else qtrend_store = NULL
    if(exo)    exo_store    = array(NA_real_, c(thindraws, bigT, n)) else exo_store = NULL
    for(irep in 1:thindraws){
      
      if(irep%%50==0) cat(paste0("HD Computation: ",irep," / ", thindraws,".\n"))
      
      temp = gen_compMat(A=run$A[irep,1:(n*plag),], n=n, p=plag)
      Cm   = temp$Cm
      Jm   = temp$Jm
      
      if(max(abs(Re(eigen(Cm)$values))) > 1)
        next
      
      Sigchol    = t(chol(run$SIGMA[irep,,]))
      SigcholR   = Sigchol %*% R
      strMat     = solve(SigcholR)
      
      # structural errors
      eps = strMat %*% t(Y - X%*%run$A[irep,,])
      
      # contribution of each shock
      invA_big = matrix(0, K, n)
      invA_big[1:n,] = strMat
      HDshock_big = array(0, c(K, bigT+1, n))
      HDshock     = array(0, c(n, bigT+1, n))
      for(nn in 1:n){
        eps_big = matrix(0, n, bigT+1)
        eps_big[nn,2:(bigT+1)] = eps[nn,]
        for(tt in 2:(bigT+1)){
          HDshock_big[,tt,nn] = invA_big %*% eps_big[,tt] + Cm %*% HDshock_big[,(tt-1),nn]
          HDshock[,tt,nn]     = t(Jm) %*% HDshock_big[,tt,nn]
        }
      }
      # [nobs x shock x var]
      hd_store[irep,,,] = aperm(HDshock[,2:(bigT+1),], c(2,3,1))
      
      # contribution of initial value
      HDinit_big = matrix(0, n*plag, bigT+1)
      HDinit     = matrix(0, n, bigT+1)
      HDinit_big[,1] = t(X[1,1:n*plag,drop=FALSE])
      HDinit[,1]     = t(Jm) %*% HDinit_big[,1]
      for(tt in 2:(bigT+1)){
        HDinit_big[,tt] = Cm %*% HDinit_big[,tt-1]
        HDinit[,tt]     = t(Jm) %*% HDinit_big[,tt]
      }
      init_store[irep,,] = t(HDinit[,2:(bigT+1)])
      
      # contribution of constant
      HDcons = matrix(0, n, bigT+1)
      if(cons){
        HDcons_big = matrix(0, n*plag, bigT+1)
        CC = matrix(0, K, 1)
        CC[1:n,] = A[irep,"cons",]
        for(tt in 2:(bigT+1)){
          HDcons_big[,tt] = CC + Cm %*% HDcons_big[,tt-1]
          HDcons[,tt]     = t(Jm) %*% HDcons_big[,tt]
        }
        cons_store[irep,,] = t(HDcons[,2:(bigT+1)])
      }
      
      # contribution of trend
      HDtrend = matrix(0, n, bigT+1)
      if(trend){
        HDtrend_big = matrix(0, K, bigT+1)
        TT = matrix(0, K, 1)
        TT[1:n,] = A[irep,"trend",]
        for(tt in 2:(bigT+1)){
          HDtrend_big[,tt] = TT*(tt-1) + Cm %*% HDtrend_big[,tt-1]
          HDtrend[,tt]     = t(Jm) %*% HDtrend_big[,tt]
        }
        trend_store[irep,,] = t(HDtrend[,2:(bigT+1)])
      }
      
      # contribution of quadratic trend
      HDqtrend = matrix(0, n, bigT+1)
      if(qtrend){
        HDqtrend_big = matrix(0, K, bigT+1)
        TT2 = matrix(0, K, 1)
        TT2[1:n,] = A[irep,"qtrend",]
        for(tt in 2:(bigT+1)){
          HDqtrend_big[,tt] = TT2*(tt-1)^2 + Cm %*% HDqtrend_big[,tt-1]
          HDqtrend[,tt]     = t(Jm) %*% HDqtrend_big[,tt]
        }
        qtrend_store[irep,,] = t(HDqtrend[,2:(bigT+1)])
      }
      
      # contributio  of exogenous regressors
      HDexo = matrix(0, n, bigT+1)
      if(exo){
        HDexo_big = matrix(0, K, bigT+1)
        EXO = matrix(0, K, Kex)
        EXO[1:n,] = A[irep,grepl("Ex",dimnames(A)),]
        for(tt in 2:(bigT+1)){
          HDexo_big[,tt] = t(EXO %*% X[,grepl("Ex",colnames(X))]) + Cm %*% HDexo_big[,tt-1]
          HDexo[,tt]     = t(Jm) %*% HDexo_big[,tt]
        }
        exo_store[irep,,] = t(HDexo[,2:(bigT+1)])
      }
    } # END-for thindraws
    
    return(list(hd_store=hd_store, init_store=init_store, cons_store=cons_store,
                trend_store=trend_store, qtrend_store=qtrend_store, exo_store=exo_store))
  }
}
