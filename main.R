###################################################################
#                                                                 #
# Replication files: Natural Gas Prices, Inflation Expectations,  #
#                    and the Pass-Through to Euro Area Inflation  #
#                                                                 #
# Maximilian Boeck                                                #
#                                                                 #
# Created:   20/09/2024                                           #
# Last Edit: 23/09/2024                                           #
###################################################################
rm(list=ls())

# CAUTION: run-time of this program can be quite extensive on a standard computer
# CAUTION: this program saves estimations on your directory (> 3 Gb)
# in case you want to delete the estimations afterward, delete the content of the following
# directory: /results/rda

# load packages
library(readxl)
library(stringr)
library(abind)
library(Hmisc)
library(stringr)
library(MASS)
library(Matrix)
library(mvtnorm)
library(tseries)

# load functions
funs <- list.files("./codes/functions_misc", full.name=TRUE)
for(fun in funs) if(grepl("R$",fun)) source(fun)

# load functions
funs <- list.files("./codes/functions_varsv", full.name=TRUE)
for(fun in funs) if(grepl("R$",fun)) source(fun)
rm(fun, funs)

# VAR settings
plag   = 12       # number of lags in the VAR
nhor   = 61       # compute impulse responses up to this horizon
hor    = 61       # impulse response horizon in plots (has to be at least <= nhor)

# MCMC settings
draws  = 25000    # number of saved draws
burnin = 10000    # number of burnin draws
thin   = 2        # thinning factor: save every thin'd draw

# Figure 1 / Figure 2
source("./codes/motivationplot.R")

# Figure 3 / Figure 4 / Figure D1a / Figure D3
setting = "baseline"
source("./codes/var_analysis.R")

# Figure 5
setting = "baseline-ils"
source("./codes/var_analysis_ILS.R")

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

# Figure D2a / Figure D2b
rob_settings = paste0("robustness-",c("core","ip","lags-six","ttf-growth"))
for(setting in rob_settings){
  source("./codes/var_analysis.R")
}
ident_settings = paste0("baseline-alternative-",c("sign1","sign2","two_shocks","one_shock"))
for(setting in ident_settings){
  source("./codes/var_analysis.R")
}
source("./codes/var_analysis_plot_robustness.R")
