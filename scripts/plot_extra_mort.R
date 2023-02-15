library(data.table)
setwd("/home/sean/Documents/zika_project/mosquito_dynamics/")
moz.dat <- fread("output/AegyptiPars_2_offset.csv")

pdf("figures/extra_mortality_ts.pdf")
plot(seq(as.Date("2000-01-01"),by="1 day",length.out=nrow(moz.dat)),
     moz.dat$extraDeath,
     bty="n",las=1,xaxs="i",yaxs="i",type='l',
     xlab="Year",ylab=expression(mu[c]))
dev.off()
