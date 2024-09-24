# specifications
vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.")
varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations ")
varAxis       = c("%", "%", "pp", "pp", "pp")
shockVar      = "Real Gas Price"
restrVar      = "Inflation Expectations "
addNames      = "Price Level"
aggNames      = "Inflation"
addAxis       = c("%")
varNames_full = c(varNames, addNames)
varAxis_full  = c(varAxis, addAxis)
n             = length(varNames)
S             = length(shockVar)
I             = length(addNames)
o             = n+I
tcode         = c(4,5,1,5,1)                         # 1=level, 2=diff, 3=2nd diff, 4=log, 5=logdiff, 6=2nd logdiff
annual        = c(0,0,0,0,0)                         # 0=not annualized, 1=annualized
scode         = c(0,0,0,1,0)                         # 0=no seas. adj; 1=seas adj.
cumul         = c(0,1,0,0,0)                         # 0=no cumulation; 1=cumulation
tcode_lag     = c(1,1,1,12,1)                        # number of lags for differences: 1=month-on-month; 12=year-on-year
freq          = 12                                   # frequency of data: 4=quarterly; 12=monthly
shock_idx     = 1                                    # position of shock variable
cf_type       = "specific"                           # 'specific': use identified idiosyncratic inflation expectations shock; 'all': use all shocks

starttime1    = c(2004,4)
starttime2    = c(2005,4)
endtime       = c(2022,12)

emp_percs     = c(.05, .10, .16, .50, .84, .90, .95)

# ils-specific
ilsNames      = c("1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "12Y", "15Y", "20Y", "25Y", "30Y")
ilsN          = length(ilsNames)

# sign-restrictions
MaxTries = 7500

# horizon
H.restr <- 1
N.restr <- H.restr*n

Smat <- Zmat <- vector(mode="list", length=5)
# first shock: real gas shock
Smat[[1]] = matrix(0, 4, N.restr)
Smat[[1]][1,1] = 1; Smat[[1]][2,3] = 1; Smat[[1]][3,4] = 1; Smat[[1]][4,5] = 1
Zmat[[1]] = matrix(0, 1, N.restr)
Zmat[[1]][1,2] <- 1
# second shock: demand shock
Smat[[2]] <- matrix(0, n, n*H.restr) 
Smat[[2]][1,1] <- 1; Smat[[2]][2,2] <- 1; Smat[[2]][3,3] <- 1; Smat[[2]][4,4] <- 1; Smat[[2]][5,5] <- 1
Zmat[[2]] <- NULL
# third shock: monetary policy shock
Smat[[3]] <- matrix(0, n, n*H.restr)
Smat[[3]][1,1] <- -1; Smat[[3]][2,2] <- -1; Smat[[3]][3,3] <- 1; Smat[[3]][4,4] <- -1; Smat[[3]][5,5] <- -1
# fourth shock: supply shock
Smat[[4]] <- matrix(0, n-1, n*H.restr)
Smat[[4]][1,1] <- -1; Smat[[4]][2,2] <- -1; Smat[[4]][3,4] <- 1; Smat[[4]][4,5] <- 1
Zmat[[4]] <- NULL
# fourth shock: inflation expectations shock
Smat[[5]] <- matrix(0, 1, n*H.restr); Smat[[5]][1,5] <- 1
Zmat[[5]] <- matrix(0, 4, n*H.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,2] <- 1; Zmat[[5]][3,3] <- 1; Zmat[[5]][4,4] <- 1

#-----------------------------------------------------------------------------------------
# Preparations
#-----------------------------------------------------------------------------------------

source("./codes/setup/00_PlotSpec.R")
source("./codes/setup/specify_var_settings.R")

# load data
ctry = "EA-m"
data = read_xlsx(path="./codes/BZ2_data.xlsx", sheet=ctry)
if(ctry == "EA-m" || ctry == "US-m") data = ts(data[,-1], start=c(2004,4), frequency=12)
if(ctry == "EA-q") data = ts(data[,-1], start=c(2004,2), frequency=4)

#-------------------------------------------------------------------------------
# container for differences
diffmax_store = array(NA_real_, c(draws, ilsN, o),       dimnames=list(NULL,       ilsNames, varNames_full))
diff_store    = array(NA_real_, c(draws, nhor, ilsN, o), dimnames=list(NULL, NULL, ilsNames, varNames_full))
maxind_store  = array(NA_real_, c(ilsN, o),              dimnames=list(            ilsNames, varNames_full))

#-------------------------------------------------------------------------------
# start big-for-loop
for(ilshor in 1:ilsN){
  # give some information
  cat(paste0("Estimation with ", ilsNames[ilshor], " Inflation Expectations.. \n"))
  
  #-----------------------------------------------------------------------------------------
  # Change Variable Names
  #-----------------------------------------------------------------------------------------
  vars.i        = vars
  varNames.i    = varNames
  
  vars.i[5]     = paste0(vars.i[5],ilsNames[ilshor])
  varNames.i[5] = paste0(varNames[5],ilsNames[ilshor])
  varNames_full = c(varNames.i, addNames)
  restrVar      = paste0("Inflation Expectations ",ilsNames[ilshor])
  
  setting.i  = paste0(setting,"-",str_extract(ilsNames[ilshor],"[0-9]+"),"y")
  filename_r = paste0("./results/rda/",ctry,"_",setting.i,"_")
  
  diffmax.i_store = array(NA_real_, c(draws, o),       dimnames=list(NULL,       varNames_full))
  diff.i_store    = array(NA_real_, c(draws, nhor, o), dimnames=list(NULL, NULL, varNames_full))
  maxind.i_store  = array(NA_real_, c(o),              dimnames=list(            varNames_full))
  #-----------------------------------------------------------------------------------------
  # Estimate Model
  #-----------------------------------------------------------------------------------------
  
  # select data
  Yraw = window(data[,vars.i], start=starttime1, end=endtime)
  Traw = nrow(Yraw)
  
  # transform data
  for(nn in 1:n){
    Yraw[,nn] = transx(Yraw[,nn], tcode[nn], lag=tcode_lag[nn])
  }
  Yraw = ts(Yraw[(max(tcode_lag)+1):nrow(Yraw),,drop=FALSE], start=starttime2, frequency=freq)
  
  # covid dummy
  Exraw = ts(matrix(0,length(Yraw),3), start=starttime2, end=endtime, frequency=freq)
  window(Exraw, start=c(2020,3), end=c(2020,3))[,1] = 1
  window(Exraw, start=c(2020,4), end=c(2020,4))[,2] = 1
  window(Exraw, start=c(2020,5), end=c(2020,5))[,3] = 1
  colnames(Exraw) = c("March 2020", "April 2020", "May 2020")
  
  # estimation
  dirName_est = paste0(filename_r,"mod_ntot=",draws+burnin,".rda")
  if(file.exists(dirName_est)){
    load(dirName_est)
  }else{
    set.seed(471)
    args = list(draws = draws, burnin = burnin, thin=thin, cons = TRUE, SV=FALSE, Ex=Exraw, kappa=kappa)
    run  = bvarsv_ng(Yraw, plag, args)
    
    # save
    save(run, file=dirName_est)
  }
  
  thindraws  = run$args$thindraws
  
  #-----------------------------------------------------------------------------------------
  # IMPULSE RESPONSE FUNCTIONS: SIGN RESTRICTIONS
  #-----------------------------------------------------------------------------------------
  # impulse response calculation and counterfactuals (wong 2015 jmcb)
  dirName_irf_sign = paste0(filename_r,"irf_sign_ntot=",draws+burnin,".rda")
  
  if(file.exists(dirName_irf_sign)){
    load(dirName_irf_sign)
  }else{
    irf_sign_store   = array(NA_real_, c(thindraws, n, n, nhor), dimnames=list(NULL, varNames.i, varNames.i, seq(0,nhor-1)))
    add_sign_store   = array(NA_real_, c(thindraws, I, n, nhor), dimnames=list(NULL, addNames, varNames.i, seq(0,nhor-1)))
    Q_store          = array(NA_real_, c(thindraws, n, n))
    counter          = rep(NA_real_, thindraws)
    for(irep in 1:thindraws){
      
      if(irep%%10==0){
        print(paste0("Round: ",irep))
        par(mfrow=c(1,1))
        if(!all(is.na(counter))) hist(counter)
      }
      
      temp    <- gen_compMat(A=run$A[irep,1:(n*plag),], n=n, p=plag)
      compMat <- temp$Cm
      Jm      <- temp$Jm
      
      if(max(abs(Re(eigen(compMat)$values))) > 1)
        next
      
      SIGMA <- run$SIGMA[irep,,]
      A0    <- t(chol(SIGMA))
      irf.restr       <- matrix(NA, N.restr, n)
      irf.restr[1:n,] <- A0
      compMati        <- compMat
      if(H.restr > 1){
        for(hh in 2:H.restr){
          irf.restr[((hh-1)*n+1):(hh*n),] <- t(Jm) %*% compMati %*% Jm %*% A0
          compMati <- compMati %*% compMat
        }
      }
      colnames(irf.restr) <- varNames
      rownames(irf.restr) <- paste(rep(varNames,H.restr),".",
                                   rep(seq(0,H.restr-1),each=n),sep="")
      
      # draw rotation matrix here
      icounter <- 0
      condall  <- 0
      no.zero.restr = all(unlist(lapply(Zmat, is.null)))
      zero.order    = order(sapply(Zmat, is.null))
      zero.order    = c(5,1,2,3,4)
      while(condall == 0 && icounter < MaxTries){
        randMat <- matrix(rnorm(n^2),n,n)
        Q <- matrix(NA_real_,n,n)
        if(no.zero.restr){
          QR <- qr(randMat)
          Q <- qr.Q(QR)
        }else{
          for(nn in 1:n){
            if(is.null(Zmat[[zero.order[nn]]])){
              Z.temp <- matrix(0, 1, N.restr)
              R      <- c()
            }else{
              Z.temp <- Zmat[[zero.order[nn]]]
              R <- Z.temp%*%irf.restr
            }
            if(nn > 1){R <- rbind(R, t(Q[,zero.order[(1:(nn-1))], drop=FALSE]))}
            
            NU  <- Null(t(R))
            x_j <- randMat[,nn,drop=FALSE]
            
            q_j <- NU%*%(t(NU)%*%x_j/sqrt(as.numeric(crossprod(t(NU)%*%x_j))))
            Q[,zero.order[nn]] <- q_j
          }
        }
        
        irf.check <- irf.restr %*% Q
        
        check_all <- 0
        for(nn in 1:n){
          S.temp <- Smat[[nn]]
          if(is.null(S.temp)) next
          nidx <- nrow(S.temp)
          check1 <- S.temp %*% sign(irf.check[,nn])
          check2 <- S.temp %*% sign(irf.check[,nn])*(-1)
          if(sum(check1)  == nidx){
            check_all <- check_all + 1
          }else if(sum(check2) == nidx){
            Q[,nn] = -Q[,nn]
            check_all <- check_all + 1
          }else{
            next
          }
        }
        if(check_all == n) condall <- 1
        icounter <- icounter + 1
      }
      counter[irep]<-icounter
      
      if(icounter == MaxTries) next
      
      shock <- A0 %*% Q
      impresp1 <- array(NA_real_, c(n, n, nhor))
      impresp1[,1:n,1] <- shock
      compMati <- compMat
      for(ihor in 2:nhor){
        impresp1[,1:n,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
        compMati <- compMati %*% compMat
      }
      
      for(nn in 1:n){
        if(cumul[nn] == 1){
          impresp1[nn,,] = t(apply(impresp1[nn,,],1,cumsum)) / ifelse(tcode_lag[nn]!=1,freq,1)
        }
        if(annual[nn] == 1){
          impresp1[nn,,] = impresp1[nn,,] * freq
        }
      }
      
      # save stuff
      irf_sign_store[irep,,,]  = impresp1
      Q_store[irep,,]          = Q
      
      # compute price/gdp level
      add_sign_store[irep,addNames,,]  = t(apply(impresp1[which(varNames %in% aggNames),,],1,function(x)cumsum(x)))/12
    }
    idx_store       = which(!is.na(irf_sign_store[,1,1,1]))
    irf_sign_store  = irf_sign_store[idx_store,,,,drop=FALSE]
    add_sign_store  = add_sign_store[idx_store,,,,drop=FALSE]
    
    irf_sign_store = abind(irf_sign_store, add_sign_store, along=2)
    irf_sign       = apply(irf_sign_store,  c(2,3,4), quantile, emp_percs, na.rm=TRUE)
    add_sign       = apply(add_sign_store,  c(2,3,4), quantile, emp_percs, na.rm=TRUE)
    irf_sign       = abind(irf_sign,  add_sign,  along=2)
    
    save(irf_sign, Q_store, file=dirName_irf_sign)
    rm(add_sign_store, add_sign)
  }
  
  #-----------------------------------------------------------------------------------------
  # Counterfactuals a la ADPRR 2021 JME
  #-----------------------------------------------------------------------------------------
  dirName_ssa_sign = paste0(filename_r,"ssa_sign_ntot=",draws+burnin,".rda")
  
  if(file.exists(dirName_ssa_sign)){
    load(dirName_ssa_sign)
  }else{
    if(cf_type == "specific"){
      irfssa_tmp         = cf_adprr(run, nhor, shockVar, restrVar, varNames.i, Q_store)
    }else if(cf_type == "all"){
      irfssa_tmp         = cf_adprr_all(run, nhor, shockVar, restrVar, varNames.i, Q_store)
    }
    irfssa_sign_store  = irfssa_tmp$irfssa_store
    
    # cumulative effects or annualization
    for(irep in 1:dim(irfssa_sign_store)[[1]]){
      for(nn in 1:n){
        if(cumul[nn] == 1){
          irfssa_sign_store[irep,nn,,] = t(apply(adrop(irfssa_sign_store[irep,nn,,,drop=FALSE],drop=c(1,2)),1,cumsum)) / ifelse(tcode_lag[nn]!=1,freq,1)
        }
        if(annual[nn] == 1){
          irfssa_sign_store[irep,nn,,] = irf_ssa_sign_store[irep,nn,,] * freq
        }
      }
    }
    
    irfssa_sign        = apply(irfssa_sign_store, c(2,3,4), quantile, c(.05, .10, .16, .50, .84, .90, .95), na.rm=TRUE)
    div_sign_store     = irfssa_tmp$div_store
    modestIntStat_sign = apply(irfssa_tmp$modInt_store, c(2,3,4), quantile, c(.05, .10, .16, .50, .84, .90, .95), na.rm=TRUE)
    
    # additional stuff
    addssa_sign_store = array(NA_real_, c(dim(irfssa_sign_store)[[1]], I, S, nhor),
                              dimnames=list(NULL, addNames, paste0(shockVar, " Shock"), seq(0,nhor-1)))
    # compute price/gdp level
    for(irep in 1:dim(irfssa_sign_store)[[1]]){
      addssa_sign_store[irep,addNames,1,] = cumsum(irfssa_sign_store[irep,aggNames,1,])/12
    }
    addssa_sign  = apply(addssa_sign_store, c(2,3,4), quantile, c(.05, .10, .16, .50, .84, .90, .95), na.rm=TRUE)
    
    irfssa_sign_store = abind(irfssa_sign_store, addssa_sign_store, along=2)
    irfssa_sign       = abind(irfssa_sign,       addssa_sign,       along=2)
    
    save(irfssa_sign, div_sign_store, modestIntStat_sign, file=dirName_ssa_sign)
    rm(addssa_sign_store, addssa_sign)
  }
  
  #-----------------------------------------------------------------------------------------
  # find maximum price level response
  #-----------------------------------------------------------------------------------------
  dirName_diff = paste0(filename_r,"diff_ntot=",draws+burnin,".rda")
  
  if(file.exists(dirName_diff)){
    load(dirName_diff)
  }else{
    # get maximum of mean response
    irf_sign_mean    = apply(irf_sign_store, c(2,3,4), mean)
    irfssa_sign_mean = apply(irfssa_sign_store, c(2,3,4), mean)
    
    for(oo in 1:o){
      maxind_ssa       = which.max(irf_sign_mean[oo,shock_idx,] - irfssa_sign_mean[oo,shock_idx,])
      for(irep in 1:length(idx_store)){
        diffmax.i_store[idx_store[irep],oo] = irf_sign_store[irep,oo,shock_idx,maxind_ssa] - irfssa_sign_store[irep,oo,shock_idx,maxind_ssa]
        for(ihor in 1:nhor){
          diff.i_store[irep,ihor,oo] = irf_sign_store[irep,oo,shock_idx,ihor] - irfssa_sign_store[irep,oo,shock_idx,ihor]
        }
      }
      maxind.i_store[oo] = maxind_ssa
    }
    
    save(diffmax.i_store, diff.i_store, maxind.i_store, file=dirName_diff)
  }
  #-----------------------------------------------------------------------------------------
  # save
  #-----------------------------------------------------------------------------------------
  diffmax_store[,ilshor,] = diffmax.i_store
  maxind_store[ilshor,]   = maxind.i_store
  diff_store[,,ilshor,]   = diff.i_store
}

diffmax_post = apply(diffmax_store, c(2,3), quantile, c(.16, .50, .86), na.rm=TRUE)
diff_post    = apply(diff_store, c(2,3,4), quantile, c(.16, .50, .86), na.rm=TRUE)

pdf("./figure5.pdf", width=10, height=6)
par(mfrow=c(1,1),mar=c(4.1,4.5,2.2,2.2))
ylim1 = range(c(0,diffmax_post[,,"Inflation"]))
plot(x = 1:ilsN, y = diffmax_post[2,,"Inflation"],
     pch = 16, type = "p", axes = FALSE, ylab="", xlab="",
     col = col.ILS, cex=1.2, ylim=ylim1)
mtext("Inflation Expectation Horizon", side=1, line=2.5, font=2, cex=1.2)
abline(h=pretty(ylim1), col="grey60")
points(x = 1:ilsN, y = diffmax_post[2,,"Inflation"], pch=16, type="p", col=col.ILS, cex=1.4)
arrows(x0=1:ilsN, y0=diffmax_post[1,,"Inflation"], 
       x1=1:ilsN, y1=diffmax_post[3,,"Inflation"], 
       col=col.ILS,
       length=0.05, angle=90, code=3, lwd=1.5)
box(lwd=2, bty="l")
axis(1, at=1:ilsN, labels=ilsNames, font=2, cex.axis=1.1, lwd=2)
axis(2, at=pretty(ylim1), labels=paste0(pretty(ylim1)," pp"),
     las=2, font=2, lwd=2, cex.axis=1.2)
dev.off()


