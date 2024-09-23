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
}else if(setting == "baseline-alternative-sign1"){
  ctry          = "EA-m"
  
  Smat[[5]]     = NULL
  Zmat[[5]]     = NULL
  Smat[[5]]     = matrix(0, 3, n*H.restr); Smat[[5]][1,2] <- 1; Smat[[5]][2,4] <- 1; Smat[[5]][3,5] <- 1
  Zmat[[5]]     = matrix(0, 2, n*H.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,3] <- 1;
}else if(setting == "baseline-alternative-sign2"){
  ctry          = "EA-m"
  
  Smat[[5]]     = NULL
  Zmat[[5]]     = NULL
  Smat[[5]]     = matrix(0, 2, n*H.restr); Smat[[5]][1,2] <- 1; Smat[[5]][2,5] <- 1
  Zmat[[5]]     = matrix(0, 3, n*H.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,3]; Zmat[[5]][3,4] <- 1;
}else if(setting == "baseline-alternative-sign3"){
  ctry          = "EA-m"
  
  # first shock: real gas shock
  Smat[[1]] = matrix(0, 5, N.restr)
  Smat[[1]][1,1] = 1; Smat[[1]][2,2] = -1; Smat[[1]][3,3] = 1; Smat[[1]][4,4] = 1; Smat[[1]][5,5] = 1
  
  Zmat <- vector(mode="list", length=5)
  Zmat[[5]] <- matrix(0, 4, N.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,2] <- 1; Zmat[[5]][3,3] <- 1; Zmat[[5]][4,4] <- 1
}else if(setting == "baseline-alternative-sign4"){
  ctry          = "EA-m"
  
  Smat[[1]] = Smat[[1]][-2,]
  Smat[[2]] = Smat[[2]][-3,]
}else if(setting == "baseline-alternative-sign5"){
  ctry          = "EA-m"
  
  Smat[[1]] = Smat[[1]]
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
}else if(setting == "extension-cons-ec"){
  ctry          = "EA-m"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "infexp_cons_ec")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "Consumer Inflation Expectations")
  varNames_full = c(varNames, addNames)
  restrVar      = "Consumer Inflation Expectations"
  
  plot_sign     = "figure6c"
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
}else if(setting == "extension-coal-1y"){
  ctry          = "EA-m"
  vars          = c("Coal_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.1Y")
  varNames      = c("Real Coal Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 1Y")
  varNames_full = c(varNames, addNames)
  shockVar      = "Real Coal Price"
  restrVar      = "Inflation Expectations 1Y"
}else if(setting == "extension-coal-5y"){
  ctry          = "EA-m"
  vars          = c("Coal_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.5Y")
  varNames      = c("Real Coal Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 5Y")
  varNames_full = c(varNames, addNames)
  shockVar      = "Real Coal Price"
  restrVar      = "Inflation Expectations 5Y"
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
}else if(setting == "rob-core"){
  ctry          = "EA-m"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CCPI", "ILS.1Y")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Core Inflation", "Inflation Expectations 1Y")
  aggNames      = "Core Inflation"
  addNames      = "Core Price Level"
  varNames_full = c(varNames, addNames)
  
  plot_sign     = "figureD1"
}else if(setting == "rob-us-core"){
  ctry          = "US-m"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CCPI", "ILS.1Y")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Core Inflation", "Inflation Expectations 1Y")
  aggNames      = "Core Inflation"
  addNames      = "Core Price Level"
  varNames_full = c(varNames, addNames)
  
  scode         = c(0,0,0,0,0)
  
  plot_sign  = "figureD3a"
}else if(setting == "rob-us-oil"){
  ctry          = "US-m"
  vars          = c("Oil_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.1Y")
  varNames      = c("Real Oil Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 1Y")
  varNames_full = c(varNames, addNames)
  shockVar      = c("Real Oil Price")
  
  scode         = c(0,0,0,0,0)
  
  plot_sign = "figureD3b"
}else if(setting == "robustness-ip"){
  ctry          = "EA-m"
  vars          = c("Gas_Worldbank_real", "IP", "SHADOWSTIR", "CPI", "ILS.1Y")
  varNames      = c("Real Gas Price", "Industrial Production", "Interest Rate", "Inflation", "Inflation Expectations 1Y")
  varNames_full = c(varNames, addNames)
}else if(setting == "robustness-stir"){
  ctry          = "EA-m"
  vars          = c("Gas_Worldbank_real", "IP", "STIR", "CPI", "ILS.1Y")
  varNames      = c("Real Gas Price", "Industrial Production", "Interest Rate", "Inflation", "Inflation Expectations 1Y")
  varNames_full = c(varNames, addNames)
}else if(setting == "robustness-lags-six"){
  ctry          = "EA-m"
  plag          = 6
}else if(setting == "robustness-ttf-growth"){
  ctry         = "EA-m"
  
  tcode         = c(5,5,1,5,1)
  cumul         = c(1,1,0,0,0)
}else if(setting == "robustness-sample1"){
  ctry          = "EA-m"
  
  starttime1    = c(2010,1)
  starttime2    = c(2011,1)
}else if(setting == "robustness-sample2"){
  ctry          = "EA-m"
  
  starttime1    = c(2012,1)
  starttime2    = c(2013,1)
}else if(setting == "extension-core"){
  ctry          = "EA-m"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.1Y", "CCPI")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 1Y", "Core Inflation")
  varAxis       = c("%", "%", "pp", "pp", "pp", "pp")
  varNames_full = c(varNames, addNames)
  varAxis_full  = c(varAxis, addAxis)
  n             = length(varNames)
  o             = n+I
  
  tcode         = c(4,5,1,5,1,5)                      
  annual        = c(0,0,0,0,0,0)
  scode         = c(0,0,0,1,0,1)
  cumul         = c(0,1,0,0,0,0)
  tcode_lag     = c(1,1,1,12,1,12)
  
  # horizon
  N.restr  = H.restr*n
  
  Smat <- Zmat <- vector(mode="list", length=n)
  # first shock: real gas shock
  Smat[[1]] = matrix(0, 4, N.restr)
  Smat[[1]][1,1] = 1; Smat[[1]][2,3] = 1; Smat[[1]][3,4] = 1; Smat[[1]][4,5] = 1
  Zmat[[1]] = matrix(0, 1, N.restr)
  Zmat[[1]][1,2] <- 1
  # second shock: demand shock
  Smat[[2]] = matrix(0, 5, N.restr) 
  Smat[[2]][1,1] = 1; Smat[[2]][2,2] <- 1; Smat[[2]][3,3] <- 1; Smat[[2]][4,4] <- 1; Smat[[2]][5,5] <- 1
  #Zmat[[2]] <- NULL
  # third shock: monetary policy shock
  Smat[[3]] <- matrix(0, 5, N.restr)
  Smat[[3]][1,1] <- -1; Smat[[3]][2,2] <- -1; Smat[[3]][3,3] <- 1; Smat[[3]][4,4] <- -1; Smat[[3]][5,5] <- -1
  # fourth shock: supply shock
  Smat[[4]] <- matrix(0, 4, N.restr)
  #Smat[[4]][1,1] <- -1; Smat[[4]][2,2] <- -1; Smat[[4]][3,3] <- 1; Smat[[4]][4,4] <- 1; Smat[[4]][5,5] <- 1 # original!
  Smat[[4]][1,1] <- -1; Smat[[4]][2,2] <- -1; Smat[[4]][3,4] <- 1; Smat[[4]][4,5] <- 1
  #Zmat[[4]] <- NULL
  # fifth shock: inflation expectations shock
  Smat[[5]] <- matrix(0, 1, N.restr); Smat[[5]][1,5] <- 1
  Zmat[[5]] <- matrix(0, 4, N.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,2] <- 1; Zmat[[5]][3,3] <- 1; Zmat[[5]][4,4] <- 1
  # sixth shock: no shock to core inflation
  Smat[[6]] <- matrix(0, 1, N.restr); Smat[[6]][1,6] <- 1
}else if(setting == "extension-energy"){
  ctry          = "EA-m"
  vars          = c("Gas_Worldbank_real", "rgdp_monthly", "SHADOWSTIR", "CPI", "ILS.1Y", "ECPI")
  varNames      = c("Real Gas Price", "Real GDP", "Interest Rate", "Inflation", "Inflation Expectations 1Y", "Energy Inflation")
  varAxis       = c("%", "%", "pp", "pp", "pp", "pp")
  varNames_full = c(varNames, addNames)
  varAxis_full  = c(varAxis, addAxis)
  n             = length(varNames)
  o             = n+I
  
  tcode         = c(4,5,1,5,1,5)                      
  annual        = c(0,0,0,0,0,0)
  scode         = c(0,0,0,1,0,1)
  cumul         = c(0,1,0,0,0,0)
  tcode_lag     = c(1,1,1,12,1,12)
  
  # horizon
  N.restr  = H.restr*n
  
  Smat <- Zmat <- vector(mode="list", length=n)
  # first shock: real gas shock
  Smat[[1]] = matrix(0, 4, N.restr)
  Smat[[1]][1,1] = 1; Smat[[1]][2,3] = 1; Smat[[1]][3,4] = 1; Smat[[1]][4,5] = 1
  Zmat[[1]] = matrix(0, 1, N.restr)
  Zmat[[1]][1,2] <- 1
  # second shock: demand shock
  Smat[[2]] = matrix(0, 5, N.restr) 
  Smat[[2]][1,1] = 1; Smat[[2]][2,2] <- 1; Smat[[2]][3,3] <- 1; Smat[[2]][4,4] <- 1; Smat[[2]][5,5] <- 1
  #Zmat[[2]] <- NULL
  # third shock: monetary policy shock
  Smat[[3]] <- matrix(0, 5, N.restr)
  Smat[[3]][1,1] <- -1; Smat[[3]][2,2] <- -1; Smat[[3]][3,3] <- 1; Smat[[3]][4,4] <- -1; Smat[[3]][5,5] <- -1
  # fourth shock: supply shock
  Smat[[4]] <- matrix(0, 4, N.restr)
  #Smat[[4]][1,1] <- -1; Smat[[4]][2,2] <- -1; Smat[[4]][3,3] <- 1; Smat[[4]][4,4] <- 1; Smat[[4]][5,5] <- 1 # original!
  Smat[[4]][1,1] <- -1; Smat[[4]][2,2] <- -1; Smat[[4]][3,4] <- 1; Smat[[4]][4,5] <- 1
  #Zmat[[4]] <- NULL
  # fifth shock: inflation expectations shock
  Smat[[5]] <- matrix(0, 1, N.restr); Smat[[5]][1,5] <- 1
  Zmat[[5]] <- matrix(0, 4, N.restr); Zmat[[5]][1,1] <- 1; Zmat[[5]][2,2] <- 1; Zmat[[5]][3,3] <- 1; Zmat[[5]][4,4] <- 1
  # sixth shock: no shock to core inflation
  Smat[[6]] <- matrix(0, 1, N.restr); Smat[[6]][1,6] <- 1
}else{
  stop("Setting not available. Please check!")
}