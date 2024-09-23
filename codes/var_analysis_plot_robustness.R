# load robustness
robspec_set = c("core","ip","lags-six","ttf-growth")
irf_sign_robspec_list = irfssa_sign_robspec_list = list()
for(rr in 1:length(robspec_set)){
  load(paste0("./results/rda/",ctry,"_robustness-",robspec_set[rr],"_irf_sign_ntot=",draws+burnin,".rda"))
  load(paste0("./results/rda/",ctry,"_robustness-",robspec_set[rr],"_ssa_sign_ntot=",draws+burnin,".rda"))
  irf_sign_robspec_list[[rr]]    = irf_sign
  irfssa_sign_robspec_list[[rr]] = irfssa_sign
}
names(irf_sign_robspec_list) = names(irfssa_sign_robspec_list) = robspec_set

robident_set = c("sign1","sign2","two_shocks","one_shock")
irf_sign_robident_list = irfssa_sign_robident_list = list()
for(rr in 1:length(robident_set)){
  load(paste0("./results/rda/",ctry,"_baseline-alternative-",robident_set[rr],"_irf_sign_ntot=",draws+burnin,".rda"))
  load(paste0("./results/rda/",ctry,"_baseline-alternative-",robident_set[rr],"_ssa_sign_ntot=",draws+burnin,".rda"))
  irf_sign_robident_list[[rr]]    = irf_sign
  irfssa_sign_robident_list[[rr]] = irfssa_sign
}
names(irf_sign_robident_list) = names(irfssa_sign_robident_list) = robident_set

# load baseline
load(paste0("./results/rda/EA-m_baseline_irf_sign_ntot=",draws+burnin,".rda"))
load(paste0("./results/rda/EA-m_baseline_ssa_sign_ntot=",draws+burnin,".rda"))

pdf(file=paste0("figureD2a.pdf"), width=12, height=3)
par(mfrow=c(1,6),mar=c(2,3.7,1.5,0.2))
for(nn in 1:o){
  ylim1=range(irf_sign[,nn,shock_idx,1:hor],
              unlist(lapply(irf_sign_robspec_list, function(l) l[4,nn,shock_idx,1:hor])),
              unlist(lapply(irfssa_sign_robspec_list, function(l) l[4,nn,shock_idx,1:hor])))
  plot.ts(irf_sign[4,nn,shock_idx,1:hor], col="black", xlab="", ylab="", ylim=ylim1, 
          main=varNames_full[nn], axes=FALSE, cex.main=ifelse(setting=="baseline",cex.title.baseline,cex.title.oth))
  abline(h=pretty(range(ylim1)), col=col.line, lwd=1)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[1,nn,shock_idx,1:hor],rev(irf_sign[7,nn,shock_idx,1:hor])),
          col = col.90, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[2,nn,shock_idx,1:hor],rev(irf_sign[6,nn,shock_idx,1:hor])),
          col = col.80, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[3,nn,shock_idx,1:hor],rev(irf_sign[5,nn,shock_idx,1:hor])),
          col = col.68, border=NA)
  lines(irf_sign[4,nn,shock_idx,1:hor], col="black", lwd=3, lty=1)
  lines(irfssa_sign[4,nn,shock_idx,1:hor], col=col.cf, lwd=3, lty=1)
  for(rr in 1:length(robspec_set)){
    lines(irf_sign_robspec_list[[rr]][4,nn,shock_idx,1:hor], col="black", lwd=2, lty=2)
    lines(irfssa_sign_robspec_list[[rr]][4,nn,shock_idx,1:hor], col=col.cf, lwd=2, lty=2)
  }
  abline(h=0, col="black", lwd=2)
  axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2)
  axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1),varAxis_full[nn]))
  box(lwd=2,bty="l")
}
dev.off()

pdf(file=paste0("figureD2b.pdf"), width=12, height=3)
par(mfrow=c(1,6),mar=c(2,3.7,1.5,0.2))
for(nn in 1:o){
  ylim1=range(irf_sign[,nn,shock_idx,1:hor],
              unlist(lapply(irf_sign_robident_list, function(l) l[4,nn,shock_idx,1:hor])),
              unlist(lapply(irfssa_sign_robident_list, function(l) l[4,nn,shock_idx,1:hor])))
  plot.ts(irf_sign[4,nn,shock_idx,1:hor], col="black", xlab="", ylab="", ylim=ylim1, 
          main=varNames_full[nn], axes=FALSE, cex.main=ifelse(setting=="baseline",cex.title.baseline,cex.title.oth))
  abline(h=pretty(range(ylim1)), col=col.line, lwd=1)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[1,nn,shock_idx,1:hor],rev(irf_sign[7,nn,shock_idx,1:hor])),
          col = col.90, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[2,nn,shock_idx,1:hor],rev(irf_sign[6,nn,shock_idx,1:hor])),
          col = col.80, border=NA)
  polygon(c(1:hor,rev(1:hor)), c(irf_sign[3,nn,shock_idx,1:hor],rev(irf_sign[5,nn,shock_idx,1:hor])),
          col = col.68, border=NA)
  lines(irf_sign[4,nn,shock_idx,1:hor], col="black", lwd=3, lty=1)
  lines(irfssa_sign[4,nn,shock_idx,1:hor], col=col.cf, lwd=3, lty=1)
  for(rr in 1:length(robident_set)){
    lines(irf_sign_robident_list[[rr]][4,nn,shock_idx,1:hor], col="black", lwd=2, lty=2)
    lines(irfssa_sign_robident_list[[rr]][4,nn,shock_idx,1:hor], col=col.cf, lwd=2, lty=2)
  }
  abline(h=0, col="black", lwd=2)
  axis(1, at=seq(1,nhor,by=12), labels=seq(0,nhor,by=12),lwd=2,font=2)
  axis(2, lwd=2, font=2, las=2, at=pretty(range(ylim1)), labels=paste0(pretty(ylim1),varAxis_full[nn]))
  box(lwd=2,bty="l")
}
dev.off()






























