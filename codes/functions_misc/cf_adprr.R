cf_adprr = function(run, nhor, shockVar, restrVar, varNames, Q_store=NULL, verbose=TRUE){
  thindraws  = run$args$thindraws
  n          = ncol(run$args$Yraw)
  S          = length(shockVar)
  plag       = run$args$plag
  
  idx_restr  = which(varNames == restrVar)
  idx_shock  = which(varNames %in% shockVar)
  
  if(idx_restr %in% idx_shock)
    stop("Cannot shock variable which we want to restrict!")
  
  if(!is.null(Q_store) && dim(Q_store)[[1]] != thindraws)
    stop("Dimension of Q_store and thindraws are not the same!")
  
  irfssa_store = array(NA_real_, c(thindraws, n, S, nhor),
                       dimnames=list(NULL, varNames, paste0(shockVar," Shock"), seq(0,nhor-1)))
  div_store    = array(NA_real_, c(thindraws, 1, S),
                       dimnames=list(NULL, c("q-divergence"), paste0(shockVar, "Shock")))
  modInt_store = array(NA_real_, c(thindraws, n, S, nhor),
                       dimnames=list(NULL, varNames, paste0(shockVar," Shock"), seq(0,nhor-1)))
  for(irep in 1:thindraws){
    
    if(irep%%50==0 && verbose) print(paste0("Round: ",irep))
    
    temp    <- gen_compMat(A=run$A[irep,1:(n*plag),], n=n, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    
    if(max(abs(Re(eigen(compMat)$values))) > 1)
      next
    
    if(is.null(Q_store)) Q <- diag(n) else Q <- Q_store[irep,,]
    
    if(any(is.na(Q))) next
    
    if(run$args$SV){
      SIGMA <- apply(run$SIGMA[irep,,,], c(2,3), median)
    }else{
      SIGMA <- run$SIGMA[irep,,]
    }
    
    # ADPRR 2021 JME
    A0inv = t(Q) %*% chol(SIGMA)
    Aplus = run$A[irep,,]%*%A0inv
    B     = vector(mode="list", length=plag)
    for(pp in 1:plag){
      B[[pp]] = run$A[irep,((pp-1)*n+1):(pp*n),]
    }
    MAT = matrix(0, n*nhor, n*nhor)
    for(ihor_row in 1:nhor){
      for(ihor_col in 1:nhor){
        idx_row <- ((ihor_row-1)*n+1):(ihor_row*n)
        idx_col <- ((ihor_col-1)*n+1):(ihor_col*n)
        if(ihor_col == ihor_row){
          MAT[idx_row,idx_col] = A0inv
        }else if(ihor_col > ihor_row){
          for(jj in 1:(ihor_col-ihor_row)){
            if(jj <= plag){
              idx_col_lag <- ((ihor_col-jj-1)*n+1):((ihor_col-jj)*n)
              MAT[idx_row,idx_col] = MAT[idx_row,idx_col] + MAT[idx_row,idx_col_lag] %*% B[[jj]]
            }
          } # END-for
        } # END-else-if
      } # END-for-ihor_col
    } # END for-ihor_row
    
    #----------------------------------------------------------#
    # 'conditional-on-observables' forecasting                 #
    #----------------------------------------------------------#
    # example restrict second variable in system
    ei             = matrix(0, n, 1)
    ei[idx_restr,] = 1
    Chat           = kronecker(diag(nhor), t(ei))
    fhat           = matrix(0, nhor, 1)
    Sigfhat        = matrix(0, nhor, nhor)
    
    #----------------------------------------------------------#
    # 'conditional-on-shocks' forecasting                      #
    #----------------------------------------------------------#
    # example: shock to first variable in the system in first period
    # example: allow second (idx_restr) shock to vary, remaining restrict to zero
    Xi1       = diag(n)[-idx_restr,]
    Xir       = kronecker(diag(nhor-1),Xi1)
    Xihat     = rbind(cbind(Xi1,matrix(0,n-1,(nhor-1)*n)),
                      cbind(matrix(0,(nhor-1)*(n-1),n),Xir))
    ghat      = matrix(0,nhor*(n-1),1)
    Sigghat   = matrix(0, nhor*(n-1), nhor*(n-1))
    
    #----------------------------------------------------------#
    # structural scenario analysis                             #
    #----------------------------------------------------------#
    for(ss in 1:S){
      # specify shock
      ghat[idx_shock[ss],1] = 1
      
      f         = rbind(fhat, ghat)
      Sig       = as.matrix(bdiag(list(Sigfhat,Sigghat)))
      C         = rbind(Chat, Xihat %*% solve(t(MAT)))
      D         = rbind(Chat%*%t(MAT), Xihat)
      
      # compute stuff
      Dinv      = try(ginv(D), silent=TRUE)
      if(is(Dinv, "try-error")) next
      mueps     = Dinv %*% f
      muy       = t(MAT) %*% mueps
      if(det(Sig)==0){
        irf_ssa = matrix(muy, nhor, n, byrow=TRUE)
      }else{
        Sigeps  = Dinv %*% Sig %*% t(Dinv) + diag(nhor*n) - Dinv %*% D %*% t(D) %*% t(Dinv)
        Sigy    = t(MAT) %*% Sigeps %*% MAT
        irf_ssa = t(rmvnorm(1, mean=as.numeric(muy), sigma=Sigy))
        irf_ssa = matrix(irf_ssa, nhor, n, byrow=TRUE)
      }
      
      kluf_shock = matrix(0,n*nhor,1)
      kluf_shock[idx_shock[ss],1] = 1
      kluf_Sigma = diag(n*nhor)
      klss_shock = mueps
      klss_Sigma = diag(n*nhor)
      
      KL = 0.5*(sum(diag(solve(klss_Sigma)%*%klss_Sigma)) + 
                  t(klss_shock - kluf_shock)%*%solve(klss_Sigma)%*%(klss_shock - kluf_shock) - 
                  n*nhor + log(det(klss_Sigma)/det(kluf_Sigma)))
      div = 0.5*(1 + sqrt(1 - exp(-(2*KL)/(n*nhor))))
      
      irfuf <- matrix(t(MAT)%*%kluf_shock, nhor, n, byrow=TRUE)
      irfss <- matrix(muy, nhor, n, byrow=TRUE)
      
      # direct effect
      dirEffect = irfuf - irfss
      # scaling
      scaling = matrix(NA_real_, nhor, n)
      for(hh in 1:nhor){
        tempsum = 0
        for(jj in 1:hh){
          tempsum = tempsum + irfuf[jj,]^2
        }
        scaling[hh,] = sqrt(tempsum)
      }
      modestIntStat = dirEffect/scaling
      
      if(div < 0.9){
        irfssa_store[irep,,ss,] = t(irf_ssa)
        div_store[irep,1,ss]    = div
        modInt_store[irep,,ss,] = t(modestIntStat)
      }
    }
  } # END-for thindraws
  idx          = which(!is.na(irfssa_store[,1,1,1]))
  thindraws    = length(idx)
  irfssa_store = irfssa_store[idx,,,,drop=FALSE]
  div_store    = div_store[idx,,,drop=FALSE]
  modInt_store = modInt_store[idx,,,,drop=FALSE]
  
  return(list(irfssa_store=irfssa_store, div_store=div_store, modInt_store=modInt_store, idx=idx))
}
