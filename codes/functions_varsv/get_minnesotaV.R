get_minnesotaV <- function(n, plag, nn, kappa, sig2, k){
  V_nn = matrix(0, k, 1)
  
  # according to BEAR Technical Guide page 12
  
  # construct V_nn
  for(jj in 1:(n*plag)){
    lag <- ceiling(jj/n)
    
    idx <- jj %% n
    idx <- ifelse(idx==0,n,idx)
    if(idx == nn){ # own lag coefficient
      V_nn[jj,1] = (kappa[1]/(lag^kappa[3]))^2
    }else{ # cross-variable lag coeefficient
      V_nn[jj,1] = (sig2[nn]/sig2[idx])*(kappa[1]*kappa[2]/(lag^kappa[3]))^2
    }
  }
  # deterministics
  if(k>(n*plag)){
    for(jj in (n*plag+1):k){
      V_nn[jj,1] = sig2[nn]*(kappa[1]*kappa[4])^2
    }
  }
  
  return(V_nn)
}
