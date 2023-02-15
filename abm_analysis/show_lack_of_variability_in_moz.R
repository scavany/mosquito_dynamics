library(viridis)
load("~/Documents/zika_project/mosquito_dynamics/abm_analysis/show_lack_of_variability_in_moz.RData",
     verbose=T)
dvec = seq(as.Date("2000-01-01"),as.Date("2010-12-31"),by="1 day")
## cols <- plasma(4, alpha=c(0.2,rep(1,2)))

tiff("../figures/moz_variability.tif",
     res=600, width=4152, height=4152/2, compression="lzw")

par(mar=c(5.1,5.1,4.1,2.1))
plot(dvec, colMeans(moz.out),
     type="l",lwd=2,xaxs="i",yaxs="i",bty="n",las=1,xlim=c(as.Date("2000-01-01"),as.Date("2011-01-01")),
     xlab="Year", ylab="",xaxt="n")
axis(1,label=seq(2000,2011,1),at=seq(as.Date("2000-01-01"),as.Date("2011-01-01"),by="1 year"))
mtext("Abundance",2,4)
polygon(c(dvec,rev(dvec)),
        c(apply(moz.out,2,function(x)quantile(x,0)),
          rev(apply(moz.out,2,function(x)quantile(x,1)))),
        border=FALSE,col=adjustcolor("red",0.2))

dev.off()
