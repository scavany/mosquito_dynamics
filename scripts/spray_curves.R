library(viridis)

t <- seq(0,150,0.01)
ulv <- ifelse(t<10,0,
       ifelse(t<11,1.5,
       ifelse(t<17,1.5*exp(-59900*(t-11)),
       ifelse(t<18,1.5,
       ifelse(t<24,1.5*exp(-59900*(t-18)),
       ifelse(t<25,1.5,
       ifelse(t<31,1.5*exp(-59900*(t-25)),0)))))))
tirs <- ifelse(t<10,0,
        ifelse(t<10+90,9,9*exp(-0.0630*(t-100))))

tiff("../figures/change_in_mortality.tif",
     res=600, width=4152, height=4152/1.5, compression="lzw")
cols <- plasma(9)
plot(t,tirs,type='l',xaxs="i",yaxs="i",las=1,bty="n",
     lwd=2,col=cols[1],
     xlab="Time (days)",ylab="Change in mortality rate (deaths/day)")
lines(t,ulv,lwd=2,col=cols[4])
lines(t,rep(0,length(t)),lwd=2,col=cols[7])
legend("topright", legend=c("TIRS","ULV", "None"), col=cols[c(1,4,7)],
       lty = c("solid"),lwd=c(2,2,2),bty="n")
dev.off()
