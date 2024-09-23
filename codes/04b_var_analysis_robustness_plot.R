###################################################################
#                                                                 #
# Main: VAR with Gas Prices Adding Monetary Policy                #
#       EA/US monthly/quarterly Data                              #
#                                                                 #
# Maximilian Boeck                                                #
#                                                                 #
# Created:   03/06/2024                                           #
# Last Edit: 04/06/2024                                           #
###################################################################
setwd("/users/mboeck/dropbox/projects/!Expectations_NaturalGas/01b_coding-revision")
rm(list=ls())

# libraries
library(readxl)
library(abind)
library(Hmisc)
library(stringr)
library(MASS)
library(Matrix)
library(mvtnorm)
library(seasonal)
library(tseries)
library(tempdisagg)

# load functions
funs <- list.files("./functions_misc", full.name=TRUE)
for(fun in funs) if(grepl("R$",fun)) source(fun)
# load functions
funs <- list.files("./functions_varsv", full.name=TRUE)
for(fun in funs) if(grepl("R$",fun)) source(fun)

# specify settings:
#           - baseline:                        EA data monthly    | D2, D3
# referee requests:
#           - baseline-alternative-sign3:      impose negative sign on real GDP growth
#           - baseline-alternative-sign4:      relax assumption that monetary authority reacts immediately to residual supply and demand shock
#           - baseline-alternative-two_shocks: only identify real natural gas price shock and idiosyncratic inflation expectations shock
#           - baseline-alternative-one_shock   only identify real natural gas price shock and use combined shocks to construct counterfactual
#           - robustness-core                  exchange inflation with core inflation
#           - robustness-ip                    exchange rgdp with IP
#           - robustness-lags-six              six lags
#           - robustness-ttf-growth            real natural gas price in growth rates
#           - robustness-sample                start after Great Financial Crisis in 2010M1

setting  = "baseline"
prior    = "ng"
# specification of kappa:
#                        kap_spec = 1: tight, kappa = [0.1 0.5 1.0 1000^2]
#                        kap_spec = 2: wide,  kappa = [0.5 1.0 1.0 1000^2]
#                        kap_spec = 3; mid,   kappa = [0.5 0.5 1.0 1000^2]
kap_spec = 2  

# general stuff
plag            = 12
draws           = 25000
burnin          = 10000
thin            = 2
SV              = FALSE
nhor            = 61
hor             = 61
scale           = 0.1
check_laglength = FALSE
save_est        = TRUE
save_plot       = TRUE

# plotting
plot_chol = plot_sign = plot_stats = plot_fe = plot_all = plot_pres = plot_ils = plot_rob = NULL
pres_sign = pres_sign_cf = NULL

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
tcode         = c(4,5,1,5,1)                         # 1=level, 2=diff, 3=2nd diff, 4=log, 5=logdiff, 6=2nd logdiff
annual        = c(0,0,0,0,0)                         # 0=not annualized, 1=annualized
scode         = c(0,0,0,1,0)                         # 0=no seas. adj; 1=seas adj.
cumul         = c(0,1,0,0,0)                         # 0=no cumulation; 1=cumulation
tcode_lag     = c(1,1,1,12,1)                        # number of lags for differences: 1=month-on-month; 12=year-on-year
freq          = 12                                   # frequency of data: 4=quarterly; 12=monthly
shock_idx     = 1                                    # position of shock variable
time_dummy    = 3                                    # add time dummeis for Covid period
cf_type       = "specific"                           # 'specific': use identified idiosyncratic inflation expectations shock; 'all': use all shocks
if(kap_spec == 1) kappa = c(.1, .5, 1, 1000^2)
if(kap_spec == 2) kappa = c(.5, 1, 1, 1000^2)
if(kap_spec == 3) kappa = c(.5, .5, 1, 1000^2)

starttime1    = c(2004,4)
starttime2    = c(2005,4)
endtime       = c(2022,12)

emp_percs     = c(.05, .10, .16, .50, .84, .90, .95)

ctry       = "EA-m"
# filename
filename_p = paste0("./results_ng/",ctry,"_",setting,"_")
filename_r = paste0("./results_ng/rda/",ctry,"_",setting,"_")
filename_o = paste0("./plots/")

source("./setup/00_PlotSpec.R")
source("./setup/specify_var_settings.R")

# load robustness
robspec_set = c("core","ip","lags-six","ttf-growth")
irf_sign_robspec_list = irfssa_sign_robspec_list = list()
for(rr in 1:length(robspec_set)){
  load(paste0("./results_ng/rda/",ctry,"_robustness-",robspec_set[rr],"_irf_sign_ntot=",draws+burnin,".rda"))
  load(paste0("./results_ng/rda/",ctry,"_robustness-",robspec_set[rr],"_ssa_sign_ntot=",draws+burnin,".rda"))
  irf_sign_robspec_list[[rr]]    = irf_sign
  irfssa_sign_robspec_list[[rr]] = irfssa_sign
}
names(irf_sign_robspec_list) = names(irfssa_sign_robspec_list) = robspec_set

robident_set = c("sign3","sign4","two_shocks","one_shock")
irf_sign_robident_list = irfssa_sign_robident_list = list()
for(rr in 1:length(robident_set)){
  load(paste0("./results_ng/rda/",ctry,"_baseline-alternative-",robident_set[rr],"_irf_sign_ntot=",draws+burnin,".rda"))
  load(paste0("./results_ng/rda/",ctry,"_baseline-alternative-",robident_set[rr],"_ssa_sign_ntot=",draws+burnin,".rda"))
  irf_sign_robident_list[[rr]]    = irf_sign
  irfssa_sign_robident_list[[rr]] = irfssa_sign
}
names(irf_sign_robident_list) = names(irfssa_sign_robident_list) = robident_set

# load baseline
load(paste0(filename_r,"irf_sign_ntot=",draws+burnin,".rda"))
load(paste0(filename_r,"ssa_sign_ntot=",draws+burnin,".rda"))

if(!is.null(plot_robspec)){
  pdf(file=paste0(filename_o,plot_robspec,".pdf"), width=12, height=3)
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
}

if(!is.null(plot_robident)){
  pdf(file=paste0(filename_o,plot_robident,".pdf"), width=12, height=3)
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
}






























