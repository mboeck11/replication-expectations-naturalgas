###################################################################
#                                                                 #
# Replication files: Natural Gas Prices, Inflation Expectations,  #
#                    and the Pass-Through to Euro Area Inflation  #
#                                                                 #
# Maximilian Boeck                                                #
#                                                                 #
# Created:   20/09/2024                                           #
# Last Edit: 21/09/2024                                           #
###################################################################
rm(list=ls())

# specify settings:
#           - baseline:                    EA data monthly    | Figure 4, 5, D1a, D4
#           - extension-spf-1y/5y:         EA data quarterly  | Figure 6a, 6b, D2b
#           - extension-oil-1y/5y:         EA data monthly    | Figure 7a, 7b
#           - extension-us-1y/5y:          US data monthly    | Figure 8a, 8b
#           - rob-us-core:                 US data monthly    | Figure D3a
#           - rob-us-oil:                  US data monthly    | Figure D3b

# extension: adding one variable at a time
#           - extension-core
#           - extension-energy

# alternative settings:
#           - baseline-alternative-sign1:      idiosyncratic inflation expectations shock with positive on impact response 
#                                              on output and inflation
#           - baseline-alternative-sign2:      idiosyncratic inflation expectations shock with positive on impact response
#                                              on output
# referee requests:
#           - baseline-alternative-sign3:      impose negative sign on real GDP growth
#           - baseline-alternative-sign4:      relax assumption that monetary authority reacts immediately to supply, demand, and natural gas price shock
#           - baseline-alternative-two_shocks: only identify real natural gas price shock and idiosyncratic inflation expectations shock
#           - baseline-alternative-one_shock   only identify real natural gas price shock and use combined shocks to construct counterfactual
#           - robustness-core                  exchange inflation with core inflation
#           - robustness-ip                    exchange rgdp with IP
#           - robustness-stir                  exchange shadow rate with STIR
#           - robustness-lags-six              six lags
#           - robustness-ttf-growth            real natural gas price in growth rates
#           - robustness-sample1               start after Great Financial Crisis in 2010M1
#           - robustness-sample2               start in 2012M1 (gas more independent from oil contracts)



# load packages
library(readxl);library(stringr);library(abind);library(Hmisc);library(stringr);library(MASS);library(Matrix);library(mvtnorm);library(tseries)

# load functions
funs <- list.files("./codes/functions_misc", full.name=TRUE)
for(fun in funs) if(grepl("R$",fun)) source(fun)

# load functions
funs <- list.files("./codes/functions_varsv", full.name=TRUE)
for(fun in funs) if(grepl("R$",fun)) source(fun)
rm(fun, funs)

# Figure 1 / Figure 2
source("./codes/motivationplot.R")

# Figure 3 / Figure 4 / Figure D1a / Figure D3
setting  = "baseline"
source("./codes/var_analysis.R")

# Figure 6a / Figure 6b / Figure D1b
setting = "extension-spf-1y"
source("./codes/var_analysis.R")
setting = "extension-spf-5y"
source("./codes/var_analysis.R")

# Figure 7a / Figure 7b
setting = "extension-oil-1y"
source("./codes/var_analysis.R")
setting = "extension-oil-5y"
source("./codes/var_analysis.R")

# Figure 8a / Figure 8b
setting = "extension-us-1y"
source("./codes/var_analysis.R")
setting = "extension-us-5y"
source("./codes/var_analysis.R")

# Figure D2