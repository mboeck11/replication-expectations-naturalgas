###################################################################
#                                                                 #
# Main: Prepare Settings for VAR analysis                         #
#                                                                 #
# Maximilian Boeck                                                #
#                                                                 #
# Created:   25/09/2023                                           #
# Last Edit: 20/05/2024                                           #
###################################################################

if(setting == "baseline"){
  ctry            = "EA-m"
  
  check_laglength = FALSE
     
  plot_sign       = "figure3"
  plot_stats      = "figure4"
  
  plot_robspec    = "figureD2a"
  plot_robident   = "figureD2b"
  plot_fe         = "figureD1a"
  plot_all        = "figureD3"
}else if(setting == "baseline-ils"){
  ctry          = "EA-m"
  
  plot_ils      = "figure5"
}else if(setting == "extension-spf-1y"){
  ctry          = "EA-q"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "rw1")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "SPF Inflation Expectations 1Y")
  varNames_full = c(varNames, addNames)
  restrVar      = "SPF Inflation Expectations 1Y"
     
  nhor          = 25
  hor           = 25
  plag          = 4
  tcode_lag     = c(1,1,1,4,1)
  freq          = 4
     
  starttime1    = c(2004,2)
  starttime2    = c(2005,2)
  endtime       = c(2022,3)
     
  plot_sign     = "figure6a"
  plot_fe       = "figureD1b"
}else if(setting == "extension-spf-5y"){
  ctry          = "EA-q"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "exp5")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "SPF Inflation Expectations 5Y")
  varNames_full = c(varNames, addNames)
  restrVar      = "SPF Inflation Expectations 5Y"
  
  nhor          = 25
  hor           = 25
  plag          = 4
  tcode_lag     = c(1,1,1,4,1)
  freq          = 4
    
  starttime1    = c(2004,2)
  starttime2    = c(2005,2)
  endtime       = c(2022,3)
    
  plot_sign     = "figure6b"
}else if(setting == "extension-oil-1y"){
  ctry          = "EA-m"
  vars          = c("Oil_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.1Y")
  varNames      = c("Real Oil Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 1Y")
  varNames_full = c(varNames, addNames)
  shockVar      = "Real Oil Price"
  restrVar      = "Inflation Expectations 1Y"
  
  plot_sign     = "figure7a"
}else if(setting == "extension-oil-5y"){
  ctry          = "EA-m"
  vars          = c("Oil_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.5Y")
  varNames      = c("Real Oil Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 5Y")
  varNames_full = c(varNames, addNames)
  shockVar      = "Real Oil Price"
  restrVar      = "Inflation Expectations 5Y"
  
  plot_sign     = "figure7b"
}else if(setting == "extension-us-1y"){
  ctry          = "US-m"
  
  scode         = c(0,0,0,0,0)
  
  plot_sign     = "figure8a"
}else if(setting == "extension-us-5y"){
  ctry          = "US-m"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.5Y")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 5Y")
  varNames_full = c(varNames, addNames)
  restrVar      = "Inflation Expectations 5Y"
  
  scode         = c(0,0,0,0,0)
  
  plot_sign     = "figure8b"
}else if(setting == "robustness-core"){
    ctry          = "EA-m"
    vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CCPI", "ILS.1Y")
    varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Core Inflation", "Inflation Expectations 1Y")
    aggNames      = "Core Inflation"
    addNames      = "Core Price Level"
    varNames_full = c(varNames, addNames)
}else if(setting == "robustness-ip"){
  ctry          = "EA-m"
  vars          = c("Gas_Worldbank_real", "IP", "SHADOWSTIR", "CPI", "ILS.1Y")
  varNames      = c("Real Gas Price", "Industrial Production", "Interest Rate", "Inflation", "Inflation Expectations 1Y")
  varNames_full = c(varNames, addNames)
}else if(setting == "robustness-lags-six"){
  ctry          = "EA-m"
  plag          = 6
}else if(setting == "robustness-ttf-growth"){
  ctry         = "EA-m"
  
  tcode         = c(5,5,1,5,1)
  cumul         = c(1,1,0,0,0)
}else if(setting == "baseline-alternative-sign1"){
  ctry          = "EA-m"
  
  # first shock: real gas shock
  Smat[[1]] = matrix(0, 5, N.restr)
  Smat[[1]][1,1] = 1; Smat[[1]][2,2] = -1; Smat[[1]][3,3] = 1; Smat[[1]][4,4] = 1; Smat[[1]][5,5] = 1
  
  Zmat <- vector(mode="list", length=5)
  Zmat[[5]] <- matrix(0, 4, N.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,2] <- 1; Zmat[[5]][3,3] <- 1; Zmat[[5]][4,4] <- 1
}else if(setting == "baseline-alternative-sign2"){
  ctry          = "EA-m"
  
  Smat[[1]] = Smat[[1]][-2,]
  Smat[[2]] = Smat[[2]][-3,]
}else if(setting == "baseline-alternative-two_shocks"){
  ctry          = "EA-m"
  
  Smat[[2]] = matrix(0, 1, N.restr); Smat[[2]][1,2] = 1
  Smat[[3]] = matrix(0, 1, N.restr); Smat[[3]][1,3] = 1
  Smat[[4]] = matrix(0, 1, N.restr); Smat[[4]][1,4] = 1
}else if(setting == "baseline-alternative-one_shock"){
  ctry          = "EA-m"
  
  cf_type       = "all"
  
  Smat <- Zmat <- vector(mode="list", length=5)
  # first shock: real gas shock
  Smat[[1]] = matrix(0, 4, N.restr)
  Smat[[1]][1,1] = 1; Smat[[1]][2,3] = 1; Smat[[1]][3,4] = 1; Smat[[1]][4,5] = 1
  Zmat[[1]] = matrix(0, 1, N.restr)
  Zmat[[1]][1,2] <- 1
  Smat[[2]] = matrix(0, 1, N.restr); Smat[[2]][1,2] = 1
  Smat[[3]] = matrix(0, 1, N.restr); Smat[[3]][1,3] = 1
  Smat[[4]] = matrix(0, 1, N.restr); Smat[[4]][1,4] = 1
  Smat[[5]] = matrix(0, 1, N.restr); Smat[[5]][1,5] = 1
}else{
  stop("Setting not available. Please check!")
}