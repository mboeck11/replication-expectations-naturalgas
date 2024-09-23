bvarsv_ng <- function(Yraw, plag = 1, args = NULL){
  #----------------------------------------INPUTS----------------------------------------------------#
  # prepare arguments
  draws = burnin = 5000; SV = FALSE; cons=TRUE; trend=FALSE; qtrend=FALSE; thin=1; Ex=NULL; save.prior=FALSE; verbose=TRUE
  if(!is.null(args)){
    for(aa in c("draws","burnin","SV","cons","trend","qtrend","thin","Ex","save.prior","verbose")){
      if(aa%in%names(args)) assign(aa, args[[aa]])
    }
  }
  arglist=list(Yraw=Yraw, plag=plag, draws=draws, burnin=burnin, SV=SV, cons=cons, trend=trend, qtrend=qtrend,
               thin=thin, Ex=Ex, save.prior=save.prior)
  
  #----------------------------------------PACKAGES--------------------------------------------------#
  require(stochvol,quietly=TRUE)
  require(GIGrvg, quietly=TRUE)
  require(Rcpp, quietly=TRUE)
  require(MASS, quietly=TRUE)
  require(mvtnorm, quietly=TRUE)
  sourceCpp("./codes/functions_varsv/do_rgig.cpp")
  
  #-------------------------------------------START--------------------------------------------------#
  Traw     = nrow(Yraw)
  n        = ncol(Yraw)
  K        = n*plag
  Ylag     = mlag(Yraw,plag)
  varNames = colnames(Yraw)
  if(is.null(varNames)) varNames = paste0("Y",seq(n))
  varNameslags = NULL
  for(pp in 1:plag) varNameslags = c(varNameslags,paste(varNames,".lag",pp,sep=""))
  colnames(Ylag) = varNameslags
  
  texo = FALSE; Mex = 0; Exraw = NULL
  if(!is.null(Ex)){
    Exraw   = Ex; Mex = ncol(Exraw)
    texo    = TRUE
    ExNames = paste0("Ex.",colnames(Exraw))
    if(is.null(ExNames)) ExNames = paste0("Ex.",seq(1,Mex))
    varNameslags = c(varNameslags, ExNames)
  }
  
  Xraw = cbind(Ylag,Exraw)
  X    = Xraw[(plag+1):nrow(Xraw),,drop=FALSE]
  Y    = Yraw[(plag+1):Traw,,drop=FALSE]
  bigT = nrow(X)
  
  if(cons){
    X                    = cbind(X,1)
    varNameslags         = c(varNameslags,"cons")
    colnames(X)[ncol(X)] = "cons"
  }
  if(trend){
    X                    = cbind(X,seq(1,bigT))
    varNameslags         = c(varNameslags,"trend")
    colnames(X)[ncol(X)] = "trend"
  }
  if(qtrend){
    X                    = cbind(X,seq(1,bigT)^2)
    varNameslags         = c(varNameslags,"qtrend")
    colnames(X)[ncol(X)] = "qtrend"
  }
  
  k = ncol(X)
  v = (n*(n-1))/2
  
  #---------------------------------------------------------------------------------------------------------
  # HYPERPARAMETERS
  #---------------------------------------------------------------------------------------------------------
  # no SV
  a_1       = 3
  b_1       = 0.03
  # SV
  Bsigma    = 1
  a0        = 5
  b0        = 1.5
  bmu       = 0
  Bmu       = 100^2
  # NG
  d_lambda    = 0.01
  e_lambda    = 0.01
  tau_start   = 0.7
  sample_tau  = TRUE
  
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv = try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv = MASS::ginv(crossprod(X))
  A_OLS  = XtXinv%*%(t(X)%*%Y)
  E_OLS  = Y - X%*%A_OLS
  S_OLS  = crossprod(E_OLS)/(bigT-k)
  
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw    = A_OLS
  SIGMA     = array(S_OLS, c(n,n,bigT))
  Em_str    = E_OLS
  L_draw    = diag(n)
  L_drawinv = solve(L_draw)
  
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # prior mean
  prmean  = 1
  a_prior = matrix(0,k,n); diag(a_prior) = prmean
  l_prior = matrix(0,n,n); l_prior[lower.tri(l_prior)] = NA_real_; diag(l_prior) = NA_real_
  
  # prior variance
  A_prior = matrix(10,k,n); A_prior[k,] = 100^2
  L_prior = matrix(10,n,n); L_prior[lower.tri(L_prior)] = NA_real_; diag(L_prior) = NA_real_
  
  # NG stuff
  lambda2_draw = matrix(0.01,      plag+1, 1, dimnames=list(c(paste0("lag.",seq(plag)),"cov"),NULL))
  tau_draw     = matrix(tau_start, plag+1, 1, dimnames=dimnames(lambda2_draw))
  tau_tuning   = matrix(.43,       plag+1, 1, dimnames=dimnames(lambda2_draw))
  tau_accept   = matrix(0,         plag+1, 1, dimnames=dimnames(lambda2_draw))
  
  # SV quantities
  Sv_draw  = matrix(-3,bigT,n)
  svdraw   = list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,bigT))
  pars_var = matrix(c(-3,.9,.2,-3),4,n,dimnames=list(c("mu","phi","sigma","latent0"),NULL))
  svl      = list()
  for(nn in 1:n) svl[[nn]] = svdraw
  
  hv        = svdraw$latent
  para      = list(mu=-3,phi=.9,sigma=.2)
  Sv_priors = specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  ntot = draws+burnin
  
  # thinning
  count      = 0
  thindraws  = draws/thin
  thin.draws = seq(burnin+1,ntot,by=thin)
  arglist    = c(arglist, thindraws=thindraws)
  
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store       = array(NA_real_, c(thindraws, k, n))
  L_store       = array(NA_real_, c(thindraws, n, n))
  res_store     = array(NA_real_, c(thindraws, bigT, n))
  # SV
  if(SV){
    Sv_store    = array(NA_real_, c(thindraws, bigT, n))
    SIGMA_store = array(NA_real_, c(thindraws, bigT, n, n))
    pars_store  = array(NA_real_, c(thindraws, 4, n))
  }else{
    Sv_store    = array(NA_real_, c(thindraws, n))
    SIGMA_store = array(NA_real_, c(thindraws, n, n))
    pars_store  = NULL
  }
  # NG
  if(save.prior){
    Aprior_store  = array(NA_real_, c(thindraws, k, n))
    Lprior_store  = array(NA_real_, c(thindraws, n, n))
    lambda2_store = array(NA_real_, c(thindraws, plag+1, 1))
    tau_store     = array(NA_real_, c(thindraws, plag+1, 1))
  }else{
    Aprior_store = Lprior_store = lambda2_store = tau_store = NULL
  }
  
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    for(nn in 1:n){
      if(nn > 1){
        X1      = cbind(-Y[,1:(nn-1),drop=FALSE], X)
        p_prior = rbind(l_prior[1:(nn-1),nn,drop=FALSE],a_prior[,nn,drop=FALSE])
        V_prior = diag(1/c(L_prior[1:(nn-1),nn],A_prior[,nn]))
      }else{
        X1      = X
        p_prior = a_prior[,nn,drop=FALSE]
        V_prior = diag(1/A_prior[,nn])
      }
      Y.nn = Y[,nn,drop=FALSE] * exp(-0.5*Sv_draw[,nn])
      X.nn = X1                * exp(-0.5*Sv_draw[,nn])
      
      V_post = try(chol2inv(chol(crossprod(X.nn) + V_prior)), silent=TRUE)
      if(is(V_post,"try-error")) V_post = try(solve(crossprod(X.nn) + V_prior), silent=TRUE)
      if(is(V_post,"try-error")) V_post = ginv(crossprod(X.nn) + V_prior)
      p_post = V_post %*% (crossprod(X.nn, Y.nn) + V_prior %*% p_prior)
      
      A_draw.nn = try(p_post + t(chol(V_post))%*%rnorm(k+nn-1), silent=TRUE)
      if(is(A_draw.nn,"try-error")) A_draw.nn = rmvnorm(1, p_post, V_post)
      if(nn>1) L_draw[1:(nn-1),nn] = A_draw.nn[1:(nn-1),]
      A_draw[,nn]  = A_draw.nn[nn:(k+nn-1),]
      Em_str[,nn]  = Y[,nn] - X1 %*% A_draw.nn
    }
    rownames(A_draw) = varNameslags
    
    #----------------------------------------------------------------------------
    # Step 2: Normal-Gamma shrinkage prior
    #----------------------------------------------------------------------------
    # Normal-Gamma for Covariances
    lambda2_draw["cov",1] = rgamma(n = 1,
                                   shape = d_lambda + tau_draw["cov",1]*v,
                                   rate  = e_lambda + 0.5*tau_draw["cov",1]*sum(L_prior[upper.tri(L_prior)]))
    #Step VI: Sample the prior scaling factors for covariances from GIG
    for(nn in 2:n){
      for(ii in 1:(nn-1)){
        temp = do_rgig1(lambda = tau_draw["cov",1] - 0.5, 
                        chi    = (L_draw[ii,nn] - l_prior[ii,nn])^2, 
                        psi    = tau_draw["cov",1]*lambda2_draw["cov",1])
        L_prior[ii,nn] = ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
      }
    }
    if(sample_tau){
      #Sample L_tau through a simple RWMH step
      tau_prop      = exp(rnorm(1,0,tau_tuning["cov",1]))*tau_draw["cov",1]
      post_tau_prop = atau_post(atau=tau_prop,          thetas=L_prior[upper.tri(L_prior)], k=v, lambda2=lambda2_draw["cov",1])
      post_tau_old  = atau_post(atau=tau_draw["cov",1], thetas=L_prior[upper.tri(L_prior)], k=v, lambda2=lambda2_draw["cov",1])
      post.diff     = post_tau_prop - post_tau_old
      post.diff     = ifelse(is.nan(post.diff),-Inf,post.diff)
      if (post.diff > log(runif(1,0,1))){
        tau_draw["cov",1]   = tau_prop
        tau_accept["cov",1] = tau_accept["cov",1] + 1
      }
      if (irep<(0.5*burnin)){
        if ((tau_accept["cov",1]/irep)>0.3)  tau_tuning["cov",1] = 1.01*tau_tuning["cov",1]
        if ((tau_accept["cov",1]/irep)<0.15) tau_tuning["cov",1] = 0.99*tau_tuning["cov",1]
      }
    }
    # Normal-Gamma for endogenous variables
    for (pp in 1:plag){
      slct.i  = grep(paste0("\\.lag",pp), rownames(A_draw))
      A.lag   = A_draw[slct.i,,drop=FALSE]
      a.prior = a_prior[slct.i,,drop=FALSE]
      A.prior = A_prior[slct.i,,drop=FALSE]
      
      n.end = nrow(A.lag)
      if(pp==1){
        lambda2_draw[pp,1] = rgamma(n     = 1,
                                    shape = d_lambda + tau_draw[pp,1]*n*n.end,
                                    rate  = e_lambda + 0.5*tau_draw[pp,1]*sum(A.prior))
      }else{
        lambda2_draw[pp,1] = rgamma(n     = 1,
                                    shape = d_lambda + tau_draw[pp,1]*n*n.end,
                                    rate  = e_lambda + 0.5*tau_draw[pp,1]*prod(lambda2_draw[1:(pp-1),1])*sum(A.prior))
      }
      for(jj in 1:n.end){
        for(nn in 1:n){
          temp = do_rgig1(lambda = tau_draw[pp,1] - 0.5,
                          chi    = (A.lag[jj,nn] - a.prior[jj,nn])^2,
                          psi    = tau_draw[pp,1]*prod(lambda2_draw[1:pp,1]))
          A.prior[jj,nn] = ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
        }
      }
      A_prior[slct.i,] <- A.prior
      if (sample_tau){
        #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
        tau_prop      = exp(rnorm(1,0,tau_tuning[pp,1]))*tau_draw[pp,1]
        post_tau_prop = atau_post(atau=tau_prop,       thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,1]), k=length(A.prior))
        post_tau_old  = atau_post(atau=tau_draw[pp,1], thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,1]), k=length(A.prior))
        post.diff     = post_tau_prop-post_tau_old
        post.diff     = ifelse(is.nan(post.diff),-Inf,post.diff)
        if (post.diff > log(runif(1,0,1))){
          tau_draw[pp,1]   = tau_prop
          tau_accept[pp,1] = tau_accept[pp,1] + 1
        }
        if (irep<(0.5*burnin)){
          if((tau_accept[pp,1]/irep)>0.3)  tau_tuning[pp,1] = 1.01*tau_tuning[pp,1]
          if((tau_accept[pp,1]/irep)<0.15) tau_tuning[pp,1] = 0.99*tau_tuning[pp,1]
        }
      }
    }
    #----------------------------------------------------------------------------
    # Step 3: Sample variances
    #----------------------------------------------------------------------------
    if(SV){
      for(nn in 1:n){
        para   = as.list(pars_var[,nn]); para$nu = Inf; para$rho=0; para$beta=0
        svdraw = svsample_fast_cpp(y=Em_str[,nn], draws=1, burnin=0, designmatrix=matrix(NA_real_),
                                   priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all",
                                   startpara=para, startlatent=Sv_draw[,nn],
                                   keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1),
                                   correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0,
                                   fast_sv=default_fast_sv)
        svl[[nn]]     = svdraw
        h_            = exp(svdraw$latent[1,])
        para$mu       = svdraw$para[1,"mu"]
        para$phi      = svdraw$para[1,"phi"]
        para$sigma    = svdraw$para[1,"sigma"]
        para$latent0  = svdraw$latent0[1,"h_0"]
        pars_var[,nn] = unlist(para[c("mu","phi","sigma","latent0")])
        Sv_draw[,nn]  = log(h_)
      }
    }else{
      for(nn in 1:n){
        S_1 = a_1 + 0.5*bigT
        S_2 = b_1 + 0.5*crossprod(Em_str[,nn])
        
        sig_eta      = 1/rgamma(1,S_1,S_2)
        Sv_draw[,nn] = log(sig_eta)
      }
    }
    
    #----------------------------------------------------------------------------
    # Step 4: store draws
    #----------------------------------------------------------------------------
    if(irep %in% thin.draws){
      count <- count+1
      
      # L-inverse
      Linv_draw = solve(L_draw)
      B_draw    = A_draw %*% Linv_draw
      
      A_store[count,,]   = B_draw
      L_store[count,,]   = Linv_draw
      res_store[count,,] = Y - X %*% B_draw
      # SV
      if(SV){
        for(tt in 1:bigT){
          SIGMA_store[count,tt,,] = t(Linv_draw)%*%diag(exp(Sv_draw[tt,]))%*%Linv_draw
        }
        Sv_store[count,,]         = Sv_draw
        pars_store[count,,]       = pars_var
      }else{
        SIGMA_store[count,,]      = t(Linv_draw)%*%diag(exp(Sv_draw[1,]))%*%Linv_draw
        Sv_store[count,]          = Sv_draw[1,]
      }
      # NG
      if(save.prior){
        Aprior_store[count,,]  = A_prior
        Lprior_store[count,,]  = L_prior
        lambda2_store[count,,] = lambda2_draw
        tau_store[count,,]     = tau_draw
      }
    }
    if(irep%%50==0 && verbose) cat("Round: ",irep, "/", ntot,"\n")
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,varNameslags,varNames)
  if(SV){
    dimnames(SIGMA_store)=list(NULL,NULL,varNames,varNames)
  }else{
    dimnames(SIGMA_store)=list(NULL,varNames,varNames)
  }
  ret <- list(Y=Y, X=X, A=A_store, SIGMA=SIGMA_store, res=res_store, 
              L=L_store, Sv=Sv_store, pars=pars_store,
              Aprior=Aprior_store, Lprior=Lprior_store, lambda2=lambda2_store, tau=tau_store,
              args=arglist)
  return(ret)
}
