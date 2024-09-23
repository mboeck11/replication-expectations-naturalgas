# general stuff
plag            = 12
draws           = 25000
burnin          = 10000
thin            = 2
SV              = FALSE
nhor            = 61
hor             = 61
scale           = 0.1
save_est        = TRUE
save_plot       = TRUE

# plotting
plot_sign = plot_stats = plot_fe = plot_all = plot_ils = NULL

# specifications
vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.1Y")
varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 1Y")
varAxis       = c("%", "%", "pp", "pp", "pp")
shockVar      = "Real Gas Price"
restrVar      = "Inflation Expectations 1Y"
addNames      = "Price Level"
aggNames      = "Inflation"
addAxis       = c("%")
varNames_full = c(varNames, addNames)
varAxis_full  = c(varAxis, addAxis)
n             = length(varNames)
S             = length(shockVar)
I             = length(addNames)
o             = n+I
tcode         = c(4,5,1,5,1)                # 1=level, 2=diff, 3=2nd diff, 4=log, 5=logdiff, 6=2nd logdiff
annual        = c(0,0,0,0,0)                # 0=not annualized, 1=annualized
cumul         = c(0,1,0,0,0)                # 0=no cumulation; 1=cumulation
tcode_lag     = c(1,1,1,12,1)               # number of lags for differences: 1=month-on-month; 12=year-on-year
freq          = 12                          # frequency of data: 4=quarterly; 12=monthly
shock_idx     = 1                           # position of shock variable
cf_type       = "specific"                  # 'specific': use identified idiosyncratic inflation expectations shock; 'all': use all shocks

starttime1    = c(2004,4)
starttime2    = c(2005,4)
endtime       = c(2022,12)

emp_percs     = c(.05, .10, .16, .50, .84, .90, .95)

# sign-restrictions
MaxTries = 7500

# horizon
H.restr  = 1
N.restr  = H.restr*n

Smat <- Zmat <- vector(mode="list", length=5)
# first shock: real gas shock
Smat[[1]] = matrix(0, 4, N.restr)
Smat[[1]][1,1] = 1; Smat[[1]][2,3] = 1; Smat[[1]][3,4] = 1; Smat[[1]][4,5] = 1
Zmat[[1]] = matrix(0, 1, N.restr)
Zmat[[1]][1,2] <- 1
# second shock: demand shock
Smat[[2]] = matrix(0, n, N.restr) 
Smat[[2]][1,1] = 1; Smat[[2]][2,2] <- 1; Smat[[2]][3,3] <- 1; Smat[[2]][4,4] <- 1; Smat[[2]][5,5] <- 1
Zmat[[2]] <- NULL
# third shock: monetary policy shock
Smat[[3]] <- matrix(0, n, N.restr)
Smat[[3]][1,1] <- -1; Smat[[3]][2,2] <- -1; Smat[[3]][3,3] <- 1; Smat[[3]][4,4] <- -1; Smat[[3]][5,5] <- -1
# fourth shock: supply shock
Smat[[4]] <- matrix(0, n-1, N.restr)
#Smat[[4]][1,1] <- -1; Smat[[4]][2,2] <- -1; Smat[[4]][3,3] <- 1; Smat[[4]][4,4] <- 1; Smat[[4]][5,5] <- 1 # original!
Smat[[4]][1,1] <- -1; Smat[[4]][2,2] <- -1; Smat[[4]][3,4] <- 1; Smat[[4]][4,5] <- 1
Zmat[[4]] <- NULL
# fifth shock: inflation expectations shock
Smat[[5]] <- matrix(0, 1, N.restr); Smat[[5]][1,5] <- 1
Zmat[[5]] <- matrix(0, 4, N.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,2] <- 1; Zmat[[5]][3,3] <- 1; Zmat[[5]][4,4] <- 1

#-----------------------------------------------------------------------------------------
# Preparations
#-----------------------------------------------------------------------------------------
source("./codes/setup/00_PlotSpec.R")
source("./codes/setup/specify_var_settings.R")

filename_p = paste0("./results/",ctry,"_",setting,"_")
filename_r = paste0("./results/rda/",ctry,"_",setting,"_")

# load data
data = read_xlsx(path="./codes/BZ2_data.xlsx", sheet=ctry)
if(ctry == "EA-m" || ctry == "US-m") data = ts(data[,-1], start=c(2004,4), frequency=12)
if(ctry == "EA-q") data = ts(data[,-1], start=c(2004,2), frequency=4)

#-----------------------------------------------------------------------------------------
# Estimate Model
#-----------------------------------------------------------------------------------------

# select data
Yraw = window(data[,vars], start=starttime1, end=endtime)
Traw = nrow(Yraw)

# transform data
for(nn in 1:n){
  Yraw[,nn] = transx(Yraw[,nn], tcode[nn], lag=tcode_lag[nn])
}
Yraw = ts(Yraw[(max(tcode_lag)+1):nrow(Yraw),,drop=FALSE], start=starttime2, frequency=freq)

# covid dummy
Exraw = ts(matrix(0,length(Yraw),3), start=starttime2, end=endtime, frequency=freq)
if(freq == 12){
  window(Exraw, start=c(2020,3), end=c(2020,3))[,1] = 1
  window(Exraw, start=c(2020,4), end=c(2020,4))[,2] = 1
  window(Exraw, start=c(2020,5), end=c(2020,5))[,3] = 1
  colnames(Exraw) = c("March 2020", "April 2020", "May 2020")
}else if(freq == 4){
  window(Exraw, start=c(2020,1), end=c(2020,1))[,1] = 1
  window(Exraw, start=c(2020,2), end=c(2020,2))[,1] = 1
  window(Exraw, start=c(2020,3), end=c(2020,3))[,1] = 1
  colnames(Exraw) = c("2020 Q1", "2020 Q2", "2020 Q3")
}

# estimation
dirName_est = paste0(str_replace(filename_r,"(?<=baseline).*","_"),"mod_ntot=",draws+burnin,".rda")
if(file.exists(dirName_est)){
  load(dirName_est)
}else{
  set.seed(471)
  args = list(draws = draws, burnin = burnin, thin=thin, cons = TRUE, SV=SV, Ex=Exraw, kappa=kappa)
  run  = bvarsv_ng(Yraw, plag, args)

  # save
  if(save_est) save(run, file=dirName_est)
}

thindraws  = run$args$thindraws

#-----------------------------------------------------------------------------------------
# IMPULSE RESPONSE FUNCTIONS
#-----------------------------------------------------------------------------------------
# impulse response calculation and counterfactuals
dirName_irf_sign = paste0(filename_r,"irf_sign_ntot=",draws+burnin,".rda")

if(file.exists(dirName_irf_sign)){
  load(dirName_irf_sign)
}else{
  irf_sign_store   = array(NA_real_, c(thindraws, n, n, nhor), dimnames=list(NULL, varNames, varNames, seq(0,nhor-1)))
  add_sign_store   = array(NA_real_, c(thindraws, I, n, nhor), dimnames=list(NULL, addNames, varNames, seq(0,nhor-1)))
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
    
    if(run$args$SV){
      SIGMA <- apply(run$SIGMA[irep,,,], c(2,3), median)
    }else{
      SIGMA <- run$SIGMA[irep,,]
    }
    
    A0 <- t(chol(SIGMA))
    irf.restr       = matrix(NA, N.restr, n)
    irf.restr[1:n,] = A0
    compMati        = compMat
    if(H.restr > 1){
      for(hh in 2:H.restr){
        irf.restr[((hh-1)*n+1):(hh*n),] = t(Jm) %*% compMati %*% Jm %*% A0
        compMati                        = compMati %*% compMat
      }
    }
    colnames(irf.restr) <- varNames
    rownames(irf.restr) <- paste(rep(varNames,H.restr),".",
                                 rep(seq(0,H.restr-1),each=n),sep="")
    
    # draw rotation matrix here
    icounter <- 0
    condall  <- 0
    no.zero.restr <- all(unlist(lapply(Zmat, is.null)))
    zero.order    <- order(sapply(Zmat, is.null))
    zero.order    <- c(5,1,2,3,4)
    if(setting == "extension-core" || setting == "extension-energy") zero.order <- c(5,1,2,3,4,6)
    if(setting == "baseline-alternative-one_shock") zero.order <- c(1,2,3,4,5)
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
            R      <- Z.temp%*%irf.restr
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
    #shock <- shock%*%diag(1/diag(shock))*scale
    impresp1         = array(NA_real_, c(n, n, nhor))
    impresp1[,1:n,1] = shock
    compMati         = compMat
    for(ihor in 2:nhor){
      impresp1[,1:n,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    
    # cumulative effects or annualization
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
    
    # compute price level
    add_sign_store[irep,addNames,,]  = t(apply(impresp1[which(varNames %in% aggNames),,],1,function(x)cumsum(x)))/freq
  }
  idx_sign        = which(!is.na(irf_sign_store[,1,1,1]))
  irf_sign_store  = irf_sign_store[idx_sign,,,,drop=FALSE]
  add_sign_store  = add_sign_store[idx_sign,,,,drop=FALSE]
  
  irf_sign  = apply(irf_sign_store,  c(2,3,4), quantile, emp_percs, na.rm=TRUE)
  add_sign  = apply(add_sign_store,  c(2,3,4), quantile, emp_percs, na.rm=TRUE)
  irf_sign  = abind(irf_sign,  add_sign,  along=2)
  
  fe_sign_store = array(NA_real_, c(dim(irf_sign_store)[1], nhor))
  for(ihor in 1:(nhor-1)){
    fe_sign_store[,ihor+1] = irf_sign_store[,aggNames,shock_idx,ihor+1] - irf_sign_store[,grep("Inflation Expectations",varNames),shock_idx,ihor]
  }
  fe_sign = apply(fe_sign_store, 2, quantile, c(.05, .10, .16, .50, .84, .90, .95), na.rm=TRUE)
  
  if(save_est) save(irf_sign, Q_store, fe_sign, file=dirName_irf_sign)
  rm(irf_sign_store, add_sign_store, add_sign)
}

# check responses
if(save_plot) pdf(file=paste0(filename_p,"irf_sign.pdf"), height=6, width=10)
par(mfrow=c(2,3),mar=c(2,3.5,1.5,0.2))
for(nn in 1:o){
  ylim1=range(irf_sign[,nn,shock_idx,1:hor])
  plot.ts(irf_sign[4,nn,shock_idx,1:hor], col="black", xlab="", ylab="", ylim=ylim1, main=varNames[nn], axes=FALSE)
  abline(h=pretty(range(ylim1)), col=col.line, lwd=1)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[1,nn,shock_idx,1:hor],rev(irf_sign[7,nn,shock_idx,1:hor])),
          col = col.90, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[2,nn,shock_idx,1:hor],rev(irf_sign[6,nn,shock_idx,1:hor])),
          col = col.80, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[3,nn,shock_idx,1:hor],rev(irf_sign[5,nn,shock_idx,1:hor])),
          col = col.68, border=NA)
  lines(irf_sign[4,nn,shock_idx,1:hor], col="black", lwd=lwd.main, lty=lty.main)
  abline(h=0, col="black", lwd=2)
  axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2)
  axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1),varAxis_full[nn]))
  box(lwd=2,bty="l")
}
if(save_plot) dev.off()

if(save_plot) pdf(file=paste0(filename_p,"irf_sign_all.pdf"), height=6, width=10)
par(mfrow=c(n,n),mar=c(3.2,3.5,2.2,2.2))
for(nn in 1:n){
  for(shock_ss in 1:5){
    ylim1=range(irf_sign[,nn,shock_ss,1:hor])
    plot.ts(irf_sign[4,nn,shock_ss,1:hor], col="black", xlab="", ylab="", ylim=ylim1, main=varNames[nn], axes=FALSE)
    abline(h=pretty(range(ylim1)), col=col.line, lwd=1)
    polygon(c(1:hor,rev(1:hor)), c(irf_sign[1,nn,shock_ss,1:hor],rev(irf_sign[7,nn,shock_ss,1:hor])),
            col = col.90, border=NA)
    polygon(c(1:hor,rev(1:hor)), c(irf_sign[2,nn,shock_ss,1:hor],rev(irf_sign[6,nn,shock_ss,1:hor])),
            col = col.80, border=NA)
    polygon(c(1:hor,rev(1:hor)), c(irf_sign[3,nn,shock_ss,1:hor],rev(irf_sign[5,nn,shock_ss,1:hor])),
            col = col.68, border=NA)
    abline(h=0, col="black", lwd=lwd.zero)
    lines(irf_sign[4,nn,shock_ss,1:hor], col="black", lwd=lwd.main, lty=lty.main)
    axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2)
    axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1),varAxis_full[nn]))
    box(lwd=2, bty="l")
  }
}
if(save_plot) dev.off()

#-----------------------------------------------------------------------------------------
# Counterfactuals (ADPRR 2021 JME)
#-----------------------------------------------------------------------------------------
dirName_ssa_sign = paste0(filename_r,"ssa_sign_ntot=",draws+burnin,".rda")
if(file.exists(dirName_ssa_sign)){
  load(dirName_ssa_sign)
}else{
  if(cf_type == "specific"){
    irfssa_tmp       = cf_adprr(run, nhor, shockVar, restrVar, varNames, Q_store)
  }else if(cf_type == "all"){
    irfssa_tmp       = cf_adprr_all(run, nhor, shockVar, restrVar, varNames, Q_store)
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
  # compute price level
  for(irep in 1:dim(irfssa_sign_store)[[1]]){
    addssa_sign_store[irep,addNames,1,] = cumsum(irfssa_sign_store[irep,aggNames,1,])/freq
  }
  addssa_sign  = apply(addssa_sign_store, c(2,3,4), quantile, c(.05, .10, .16, .50, .84, .90, .95), na.rm=TRUE)
  
  irfssa_sign = abind(irfssa_sign, addssa_sign, along=2)
  
  if(save_est) save(irfssa_sign, div_sign_store, modestIntStat_sign, file=dirName_ssa_sign)
  rm(irfssa_sign_store, addssa_sign_store, addssa_sign)
}

# check responses
if(save_plot) pdf(file=paste0(filename_p,"irf_sign_ssa.pdf"), height=6, width=10)
par(mfrow=c(2,3),mar=c(2,3.5,1.5,0.2))
for(nn in 1:o){
  ylim1=range(irf_sign[,nn,shock_idx,1:hor])
  plot.ts(irf_sign[4,nn,shock_idx,1:hor], col="black", xlab="", ylab="", ylim=ylim1, main=varNames_full[nn], axes=FALSE)
  abline(h=pretty(range(ylim1)), col=col.line, lwd=1)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[1,nn,shock_idx,1:hor],rev(irf_sign[7,nn,shock_idx,1:hor])),
          col = col.90, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[2,nn,shock_idx,1:hor],rev(irf_sign[6,nn,shock_idx,1:hor])),
          col = col.80, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[3,nn,shock_idx,1:hor],rev(irf_sign[5,nn,shock_idx,1:hor])),
          col = col.68, border=NA)
  lines(irf_sign[4,nn,shock_idx,1:hor], col="black", lwd=lwd.main, lty=lty.main)
  lines(irfssa_sign[4,nn,shock_idx,1:hor], col=col.cf, lwd=lwd.cf, lty=lty.cf)
  abline(h=0, col="black", lwd=2)
  axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2)
  axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1),varAxis_full[nn]))
  box(lwd=2,bty="l")
}
if(save_plot) dev.off()

#-----------------------------------------------------------------------------------------
# PLOTS for PAPER
#-----------------------------------------------------------------------------------------

if(!is.null(plot_sign)){
  if(setting == "baseline"){
    pdf(file=paste0(plot_sign,".pdf"), width=12, height=6)
    par(mfrow=c(2,3),mar=c(2,3.7,1.5,0.2))
  }else{
    pdf(file=paste0(plot_sign,".pdf"), width=12, height=3)
    par(mfrow=c(1,6),mar=c(2,3.7,1.5,0.2))
  }
  for(nn in 1:o){
    ylim1=range(irf_sign[,nn,shock_idx,1:hor])
    plot.ts(irf_sign[4,nn,shock_idx,1:hor], col="black", xlab="", ylab="", ylim=ylim1, 
            main=varNames_full[nn], axes=FALSE, cex.main=ifelse(setting=="baseline",cex.title.baseline,cex.title.oth))
    abline(h=pretty(range(ylim1)), col=col.line, lwd=1)
    polygon(c(1:hor,rev(1:hor)), c(irf_sign[1,nn,shock_idx,1:hor],rev(irf_sign[7,nn,shock_idx,1:hor])),
            col = col.90, border=NA)
    polygon(c(1:hor,rev(1:hor)), c(irf_sign[2,nn,shock_idx,1:hor],rev(irf_sign[6,nn,shock_idx,1:hor])),
            col = col.80, border=NA)
    polygon(c(1:hor,rev(1:hor)), c(irf_sign[3,nn,shock_idx,1:hor],rev(irf_sign[5,nn,shock_idx,1:hor])),
            col = col.68, border=NA)
    lines(irf_sign[4,nn,shock_idx,1:hor], col="black", lwd=lwd.main, lty=lty.main)
    lines(irfssa_sign[4,nn,shock_idx,1:hor], col=col.cf, lwd=lwd.cf, lty=lty.cf)
    abline(h=0, col="black", lwd=2)
    axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2)
    axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1),varAxis_full[nn]))
    box(lwd=2,bty="l")
  }
  dev.off()
}

if(!is.null(plot_stats)){
  pdf(file=paste0(plot_stats,".pdf"), width=12, height=4)
  par(mfrow=c(1,2), mar=c(2,3.5,2.2,2.2))
  ylim1 = c(-2,2)
  plot.ts(modestIntStat_sign[4,restrVar,1,], axes=FALSE, xlab="", ylab="", col="black", ylim=ylim1,
          main="SSA with Inflation Expectations Shocks", cex.main=1.5)
  abline(h=seq(-2,2,by=1), col="grey55", lwd=2, lty=2)
  polygon(c(1:hor,rev(1:hor)), c(modestIntStat_sign[1,restrVar,1,],rev(modestIntStat_sign[7,restrVar,1,])),
          col = col.90, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(modestIntStat_sign[2,restrVar,1,],rev(modestIntStat_sign[6,restrVar,1,])),
          col = col.80, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(modestIntStat_sign[3,restrVar,1,],rev(modestIntStat_sign[5,restrVar,1,])),
          col = col.68, border=NA)
  abline(h=0, col="black", lwd=lwd.zero)
  lines(modestIntStat_sign[4,restrVar,1,], col="black", lwd=lwd.main, lty=lty.main)
  axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2,cex.axis=cex.axis.4var)
  axis(2, at=seq(-2,2,by=1), lwd=2, font=2, las=2, cex.axis=cex.axis.4var, las=2)
  box(lwd=2)
  
  hist(div_sign_store, freq=FALSE, breaks=100, xlab="", ylab="", col="grey45", axes=FALSE,
       main="q-Divergence", cex.main=1.5)
  axis(1, lwd=2, font=2, cex.axis=cex.axis.4var)
  axis(2, lwd=2, font=2, las=2, cex.axis=cex.axis.4var)
  box(lwd=2,bty="l")
  abline(v=mean(div_sign_store), col="black", lwd=2)
  text(x=0.523, y=80, paste0("mean = ", round(mean(div_sign_store),2)), font=2)
  box(lwd=2)
  dev.off()
}

if(!is.null(plot_all)){
  pdf(file=paste0(plot_all,".pdf"), width=20, height=12)
  par(mar=c(2,4.5,2,2))
  par(fig=c(0.00,0.20,0.90,1.00))
  plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1), xlab="", ylab="")
  text(0.50, 0.50, "Real Natural Gas Price Shock", cex=1.1, font=2)
  par(fig=c(0.20,0.40,0.90,1.00),new=TRUE)
  plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1), xlab="", ylab="")
  text(0.50, 0.50, "Demand Shock", cex=1.1, font=2)
  par(fig=c(0.40,0.60,0.90,1.00),new=TRUE)
  plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1), xlab="", ylab="")
  text(0.50, 0.50, "Monetary Policy Shock", cex=1.1, font=2)
  par(fig=c(0.60,0.80,0.90,1.00),new=TRUE)
  plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1), xlab="", ylab="")
  text(0.50, 0.50, "Residual Supply Shock", cex=1.1, font=2)
  par(fig=c(0.80,1.00,0.90,1.00),new=TRUE)
  plot(-10,-10,axes=FALSE,ylim=c(0,1),xlim=c(0,1), xlab="", ylab="")
  text(0.50, 0.50, "Idiosyncratic Inf. Exp Shock", cex=1.1, font=2)
  share <- floor(0.90/n*1000)/1000
  for(nn in 1:n){
    for(shock_ss in 1:n){
      par(fig=c(0.2*(shock_ss-1),0.2*shock_ss,0.93-share*nn,0.93-share*(nn-1)),new=TRUE)
      ylim1=range(irf_sign[,nn,shock_ss,1:hor])
      plot.ts(irf_sign[4,nn,shock_ss,1:hor], col="black", xlab="", ylab="", ylim=ylim1, main=varNames[nn], axes=FALSE)
      abline(h=pretty(range(ylim1)), col=col.line, lwd=0.8)
      abline(v=seq(1,nhor,by=12), col=col.line, lwd=0.8)
      polygon(c(1:hor,rev(1:hor)), c(irf_sign[1,nn,shock_ss,1:hor],rev(irf_sign[7,nn,shock_ss,1:hor])), col = col.90, border=NA)
      polygon(c(1:hor,rev(1:hor)), c(irf_sign[2,nn,shock_ss,1:hor],rev(irf_sign[6,nn,shock_ss,1:hor])), col = col.80, border=NA)
      polygon(c(1:hor,rev(1:hor)), c(irf_sign[3,nn,shock_ss,1:hor],rev(irf_sign[5,nn,shock_ss,1:hor])), col = col.68, border=NA)
      abline(h=0, col="black", lwd=lwd.zero)
      lines(irf_sign[4,nn,shock_ss,1:hor], col="black", lwd=lwd.main, lty=lty.main)
      axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2)
      axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1),varAxis_full[nn]))
      box(lwd=2, bty="l")
      }
  }
  dev.off()
}

if(!is.null(plot_fe)){
  pdf(file=paste0(plot_fe,".pdf"), width=8, height=6)
  par(mfrow=c(1,1), mar=c(3,5,2,2))
  ylim1=range(fe_sign,na.rm=TRUE)
  plot.ts(fe_sign[4,1:hor], col="black", xlab="", ylab="", 
          ylim=ylim1, main="1Y Ahead Inflation Forecast Errors", axes=FALSE, cex.main=1.7)
  abline(h=pretty(range(ylim1)), col=col.line, lwd=1)
  polygon(c(1:hor,rev(1:hor)), c(fe_sign[1,1:hor],rev(fe_sign[7,1:hor])),
          col = col.90, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(fe_sign[2,1:hor],rev(fe_sign[6,1:hor])),
          col = col.80, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(fe_sign[3,1:hor],rev(fe_sign[5,1:hor])),
          col = col.68, border=NA)
  abline(h=0, col="black", lwd=lwd.zero)
  lines(fe_sign[4,1:hor], col="black", lwd=lwd.main, lty=lty.main)
  axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2,cex.axis=cex.axis.4var)
  axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1), "pp"), cex.axis=cex.axis.4var)
  box(lwd=2,bty="l")
  dev.off()
}

