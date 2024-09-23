###################################################################
#                                                                 #
# Main: Prepare Data for VAR analysis                             #
#                                                                 #
# Maximilian Boeck                                                #
#                                                                 #
# Created:   24/09/2023                                           #
# Last Edit: 24/09/2023                                           #
###################################################################
rm(list=ls())
#setwd("/users/mboeck/dropbox/projects/")
setwd("/Users/piaheckl/Dropbox/!Expectations and Natural Gas/replication-files")

# load data
library(xlsx)
library(stringr)
library(tempdisagg)
library(readxl)
library(seasonal)

#-----------------------------------------------------------------#
# data for EA-m                                                   #
#-----------------------------------------------------------------#
# macro data
data = read_xlsx(path = "../02_data/Data_monthly.xlsx", sheet = "EA")
data$TTF_IMF_real       = data$TTF_IMF*(1/data$FXUSEUR)/data$HICP
data$Gas_IMF_real       = data$Gas_IMF1*(1/data$FXUSEUR)/data$HICP
data$Gas_Worldbank_real = data$Gas_Worldbank*(1/data$FXUSEUR)/data$HICP
data$Oil_real           = data$DCOILBRENTEU/data$HICP
data$Coal_real          = data$Coal*(1/data$FXUSEUR)/data$HICP
data = ts(data[,2:ncol(data)], start=c(1994,1), frequency=12)
data = window(data, start=c(2004,4), end=c(2022,12))

# deasonalize
data[,"HICP"]        = final(seas(data[,"HICP"]))
data[,"COREHICP"]    = final(seas(data[,"COREHICP"]))
data[,"HICP_energy"] = final(seas(data[,"HICP_energy"]))

# standardise names
colnames(data)[colnames(data) == "HICP"]        = "CPI"
colnames(data)[colnames(data) == "COREHICP"]    = "CCPI"
colnames(data)[colnames(data) == "HICP_energy"] = "ECPI"

# load quarterly data
data_quarterly = read_xlsx(path="../02_data/Data_quarterly.xlsx", sheet="EA")
data_quarterly = ts(data_quarterly[,2:ncol(data_quarterly)], start=c(1999,1), frequency=4)
data_quarterly = window(data_quarterly, start=c(2004,2), end=c(2022,4))
temp_disagg    = td(data_quarterly[,"RGDP"] ~ data[,"IP"])
rgdp_monthly   = predict(temp_disagg)

# ILS data
load("../02_data/ILS_monthly.rda")
ILS = ts(as.matrix(ILS_monthly[["EA"]]), start=c(2004,4), frequency=12)
ILS = window(ILS, start=c(2004,4), end=c(2022,12))
colnames(ILS)<-paste0("ILS.",colnames(ILS))

data = ts.union(data, rgdp_monthly, ILS)
colnames(data) = str_remove(colnames(data),"(data|ILS)\\.")
rm(data_quarterly, temp_disagg, ILS_monthly, ILS_NS_monthly, rgdp_monthly, ILS)

rownames(data) = as.character(seq.Date(as.Date("2004-04-01"),as.Date("2022-12-01"),by="1 month"))

write.xlsx(data, file="./codes/BZ2_data.xlsx", sheetName="EA-m", row.names=TRUE, col.names=TRUE)
rm(data)

#-----------------------------------------------------------------#
# data for EA-q                                                   #
#-----------------------------------------------------------------#

# macro data
data = read_xlsx(path = "../02_data/Data_monthly.xlsx", sheet = "EA")
data$TTF_IMF_real       = data$TTF_IMF*(1/data$FXUSEUR)/data$HICP
data$Gas_IMF_real       = data$Gas_IMF1*(1/data$FXUSEUR)/data$HICP
data$Gas_Worldbank_real = data$Gas_Worldbank*(1/data$FXUSEUR)/data$HICP
data$Oil_real           = data$DCOILBRENTEU/data$HICP
data$Coal_real          = data$Coal*(1/data$FXUSEUR)/data$HICP
data                    = ts(data[,2:ncol(data)], start=c(1994,1), frequency=12)
data                    = window(data, start=c(2004,4), end=c(2022,12))
data                    = aggregate(data, nfrequency=4, FUN=mean)

# deasonalize
data[,"HICP"]        = final(seas(data[,"HICP"]))
data[,"COREHICP"]    = final(seas(data[,"COREHICP"]))
data[,"HICP_energy"] = final(seas(data[,"HICP_energy"]))

# standardise names
colnames(data)[colnames(data) == "HICP"]        = "CPI"
colnames(data)[colnames(data) == "COREHICP"]    = "CCPI"
colnames(data)[colnames(data) == "HICP_energy"] = "ECPI"

# macro quarterly data
data_quarterly           = read_xlsx(path="../02_data/Data_quarterly.xlsx", sheet="EA")
data_quarterly           = ts(data_quarterly[,2:ncol(data_quarterly)], start=c(1999,1), frequency=4)
data_quarterly           = window(data_quarterly[,c("RGDP","EA_Infexp_long")], start=c(2004,2))
colnames(data_quarterly) = c("rgdp_monthly","EA_infexp_long")

# ILS data
load("../02_data/ILS_quarterly.rda")
ILS = ts(as.matrix(ILS_quarterly[["EA"]]), start=c(2004,2), frequency=4)
ILS = window(ILS, start=c(2004,2), end=c(2022,4))
colnames(ILS)<-paste0("ILS.",colnames(ILS))

# survey of professional forecasters
load("../02_data/ECB SPF/SPF_ECB_ts.rda")
spf_inf <- window(spf_inf[,c("exp0","rw1","exp1","exp5","sd_exp0","sd_rw1","sd_exp1","sd_exp5")], c(2004,2), end=c(2022,4))
rm(spf_grw, spf_une)
spf_inf[,"rw1"]    <- c(spf_inf[-1,"rw1"],NA)    # rw1 is leading
spf_inf[,"sd_rw1"] <- c(spf_inf[-1,"sd_rw1"],NA)

data = ts.union(data, data_quarterly, ILS, spf_inf)
colnames(data) = str_remove(colnames(data),"(data|data_quarterly|ILS|spf_inf)\\.")
rm(ILS_quarterly, ILS, ILS_NS_quarterly, data_quarterly, spf_inf)

rownames(data) = as.character(seq.Date(as.Date("2004-04-01"),as.Date("2022-12-01"),by="1 quarter"))
write.xlsx(data, file="./codes/BZ2_data.xlsx", sheetName="EA-q", append=TRUE, row.names=TRUE, col.names=TRUE)
rm(data)

#-----------------------------------------------------------------#
# data for US-m                                                   #
#-----------------------------------------------------------------#

# macro data
data = read_xlsx(path = "../02_data/Data_monthly.xlsx", sheet = "US")
data$Oil_real           = data$WTISPLC/data$CPIAUCSL
data$Gas_IMF_real       = data$Gas_IMF2/data$CPIAUCSL
data$Gas_Worldbank_real = data$Gas_Worldbank/data$CPIAUCSL
data$Coal_real          = data$Coal/data$CPIAUCSL
data                    = ts(data[,2:ncol(data)], start=c(1990,1), frequency=12)
data                    = window(data, start=c(2004,4), end=c(2022,12))

# standardize names
colnames(data)[colnames(data) == "CPIAUCSL"] = "CPI"
colnames(data)[colnames(data) == "INDPRO"]   = "IP"
colnames(data)[colnames(data) == "FEDFUNDS"] = "STIR"
colnames(data)[colnames(data) == "CORECPI"]  = "CCPI"

# macro quarterly data
data_quarterly = read_xlsx(path="../02_data/Data_quarterly.xlsx", sheet="US")
data_quarterly = ts(data_quarterly[,2:ncol(data_quarterly)], start=c(1947,1), frequency=4)
data_quarterly = window(data_quarterly, start=c(2004,2), end=c(2022,4))
temp_disagg    = td(data_quarterly[,"GDPC1"] ~ data[,"IP"])
rgdp_monthly   = predict(temp_disagg)

# ILS data
load("../02_data/ILS_monthly.rda")
ILS = ts(as.matrix(ILS_monthly[["US"]]), start=c(2004,4), frequency=12)
ILS = window(ILS, start=c(2004,4), end=c(2022,12))
colnames(ILS)<-paste0("ILS.",colnames(ILS))

data = ts.union(data, rgdp_monthly, ILS)
colnames(data) = str_remove(colnames(data),"(data|data_quarterly|ILS)\\.")
rm(ILS_monthly, data_quarterly, ILS_NS_monthly, ILS, temp_disagg, rgdp_monthly)

rownames(data) = as.character(seq.Date(as.Date("2004-04-01"),as.Date("2022-12-01"),by="1 month"))
write.xlsx(data, file="./codes/BZ2_data.xlsx", sheetName="US-m", append=TRUE, row.names=TRUE, col.names=TRUE)
rm(data)

