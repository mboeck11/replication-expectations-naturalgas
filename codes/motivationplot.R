# load data
data_ea = read_xlsx(path="./codes/BZ2_data.xlsx", sheet="EA-m")
data_ea = ts(data_ea, start=c(2004,4), frequency=12)

#-----------------------------------------------------------------------------------------#
# FIGURE 1: Real Gas Prices and Inflation Expectations in the Euro Area.                  #
#-----------------------------------------------------------------------------------------#

# paper
pdf(file = "./figure1.pdf", width=10, height=5)
par(mar=c(3,3,3,3) + 0.5, mfrow=c(1,1))
plot(data_ea[,"Gas_Worldbank_real"], col="black", lwd=3, axes=FALSE, xlab="", ylab="")
axis(side = 2, at=pretty(range(data_ea[,"Gas_Worldbank_real"],na.rm=TRUE)), col="black", col.axis="black",las=1,font=2,lwd=3)
text(x=2005.5,y=0.47, "EA Real Gas Price (EUR/mmbtu)", xpd=NA,font=2)
box(lwd=3)
abline(h=pretty(range(data_ea[,"Gas_Worldbank_real"],na.rm=TRUE)), col="grey60")
abline(v=seq(2004,2022,by=2), col="grey60")
lines(data_ea[,"Gas_Worldbank_real"], col="black", lwd=3)
par(new=TRUE)
plot(data_ea[,"ILS.1Y"], col="#ab2830", lwd=3, axes=FALSE, xlab="", ylab="",
     ylim=range(data_ea[,"ILS.1Y"]), lty=6)
lines(data_ea[,"ILS.10Y"], col="#425df5", lty=5, lwd=3)
axis(side = 4, at=pretty(range(data_ea[,c("ILS.1Y","ILS.10Y")],na.rm=TRUE)), col="black",col.axis="black",las=1,font=2,lwd=3)
text(x=2022,y=7.9, "Inflation Linked Swaps", xpd=NA,font=2,col="black")
axis(side = 1, at=seq(2004,2022,by=2),lwd=3,col.axis="black",las=1,font=2)
legend("topleft", c("Real Gas Price", "Inflation Linked Swaps 1Y", "Inflation Linked Swaps 10Y"), col=c("black","#ab2830", "#425df5"),
       lwd=3, lty=c(1,6,5), bty="n")
dev.off()

#-----------------------------------------------------------------------------------------#
# FIGURE 2: Real Natural Gas Prices and Real Oil Prices.                                  #                                                      #
#-----------------------------------------------------------------------------------------#

commodities_ea <- data_ea[,c("Gas_Worldbank","DCOILBRENTEU")]
commodities_ea <- ts(apply(commodities_ea,2,scale),start=c(2004,1),frequency=12)

# paper
pdf(file = "./figure2.pdf", width=10, height=5)
par(mar=c(3,3,3,3) + 0.3, mfrow=c(1,1))
plot.ts(commodities_ea[,"DCOILBRENTEU"], col="#ab2830", lwd=3, axes=FALSE, xlab="", ylab="",
        ylim=range(commodities_ea), lty=5)
axis(side = 2, at=pretty(range(commodities_ea)), col="black", col.axis="black",las=1,font=2,lwd=2,cex=1.1)
box(lwd=3)
abline(h=pretty(range(commodities_ea)), col="grey60")
abline(v=seq(2004,2022,by=2), col="grey60")
lines(commodities_ea[,"DCOILBRENTEU"], col="#ab2830", lwd=3, lty=5)
lines(commodities_ea[,"Gas_Worldbank"],col="black", lwd=3)
axis(side = 1, at=seq(2004,2022,by=2),lwd=2,col.axis="black",las=1,font=2,cex=1.1)
legend("topleft", c("Natural Gas","Oil"), col=c("black","#ab2830"), lwd=3, lty=c(1,5), bty="n", cex=1.4)
dev.off()

rm(commodities_ea)
rm(data_ea)

