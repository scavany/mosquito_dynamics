setwd("~/Documents/zika_project/outbreak_response/IquitoSim_code/output")
mozzes <- read.csv("sim_nospray_immature.csv")
mozzes.gam <- read.csv("../../../mosquito_dynamics/output/mosquito_ts.csv")
dates <- 2000 + (seq(length(mozzes[,"Adults"]))/365)
cols <- c(red="#AE1A19", black="#252525", grey="#B1B1A6",
          brown="#604B4A", lgrey="#D0D0D0", white="#FFFFFF")
n.locs <- dim(read.csv("../input/Locations20190109.csv"))[1]

pdf("mosquito_abundance.pdf")
plot(dates,
     mozzes.gam[,1],
     ylab = "Total mosquito abundance",
     xlab = "Year",
     ylim=c(0,854541),
     col=cols["grey"], , cex.axis=1.3, cex.lab=1.3)
lines(dates, mozzes[,"Adults"], col=cols['red'])
legend(2006, 8e5, , legend=c("GAM", "ABM, no mvmt"),
       col=cols[c("grey", "black")], pch=c("o", NA),
       lty=c(NA, "solid"), cex=1.3)
dev.off()

pdf("zero_locations.pdf")
plot(dates,
     mozzes[,'Zero.Locations']/n.locs,
     type='l',
     ylab = "Proportion locations without mosquitoes",
     xlab = "Year",
     col=cols["black"], cex.axis=1.3, cex.lab=1.3)
dev.off()


mozzes.new <- mozzes[,2:5]
mozzes.totals <- apply(mozzes.new, 1, sum)
moz.e <- mozzes.new[,1]/mozzes.totals
moz.el <- apply(mozzes.new[,1:2],1, sum)/mozzes.totals
moz.elp <- apply(mozzes.new[,1:3],1, sum)/mozzes.totals

pdf("population_types.pdf", width=14, height=7)
plot(x=dates, y=rep(-1, length(dates)),
     col=cols["grey"], type='l', ylim=c(0,1),
     xlab="Year", ylab="Proportion", cex.axis=1.3, cex.lab=1.3)
polygon(x=c(dates, rev(dates)),
        y=c(rep(0, length(dates)), rev(moz.e)),
        border=NA, col=cols["red"])
polygon(x=c(dates, rev(dates)),
        y=c(moz.e, rev(moz.el)),
        border=NA, col=cols["grey"])
polygon(x=c(dates, rev(dates)),
        y=c(moz.el, rev(moz.elp)),
        border=NA, col=cols["brown"])
polygon(x=c(dates, rev(dates)),
        y=c(moz.elp, rep(1, length(dates))),
        border=NA, col=cols["lgrey"])
legend(x=2008, y=0.4, legend=rev(c("Eggs", "Larvae", "Pupae", "Adults")),
       col=rev(cols[c("red", "grey", "brown", "lgrey")]),
       fill=rev(cols[c("red", "grey", "brown", "lgrey")]),cex=1.3)
dev.off()
