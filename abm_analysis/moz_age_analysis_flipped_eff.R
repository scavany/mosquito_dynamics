library(fields)
setwd("~/Documents/zika_project/mosquito_dynamics/abm_analysis")
sim.nums <- c(tirs="00000", ulv="00010", none="00020")
sim.names <- names(sim.nums)
spatial.patterns <- array(dim=c(3, 41, 35)) 
dates <- seq(as.Date("2000-01-01"), as.Date("2010-12-31"), "1 day")[seq(1,4018,100)]
for (sim.num in sim.nums) {
    k <- which(sim.nums==sim.num)
    for (i in seq(41)) {
        day.num <- toString((i-1)*100)
        temp.data <- read.csv(paste0("moz_age_dist_het_flipped_eff/heterogeneity_sim_",
                                     sim.num, "_day_", day.num, ".csv"))
        temp.data$zoneID <- as.numeric(substring(temp.data$locID, 1, 2))
        for (j in seq(35)) {
            spatial.patterns[k,i,j] <- sum(temp.data$mozzes[which(temp.data$zoneID==j)])
        }
    }
}


for (sim.num in sim.nums) {
    k <- which(sim.nums==sim.num)
    spatial.patterns.sum <- apply(spatial.patterns[k,,], 1, sum)
    spatial.patterns.normalized <- spatial.patterns[k,,]/spatial.patterns.sum
    tiff(paste0("figures/spatial_patterns_normalized_flipped_eff_",
                names(sim.nums)[k],".tif"),
         res=600, width=4152, height=4152, compression="lzw")
    if(which(sim.nums==sim.num) != "none"){
        if(which(sim.nums==sim.num) == "ulv"){
            image.plot(x=dates, y=seq(1,35), z=spatial.patterns.normalized,
                       col=plasma(200), ylab="Zone", xlab="Year",zlim=c(0,max(spatial.patterns.normalized)),
                       main="High efficacy, low residuality")
        } else {
            image.plot(x=dates, y=seq(1,35), z=spatial.patterns.normalized,
                       col=plasma(200), ylab="Zone", xlab="Year",zlim=c(0,max(spatial.patterns.normalized)),
                       main="Low efficacy, high residuality")
        }
        abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
    }
    dev.off()
}

library(imager)
tiff("./figures/combined_flipped_eff.tif",compression="lzw",res=600,unit="in",width=7,height=3.5)
par(mfrow=c(1,2))
par(mar = c(0,0,0,0))
im <- load.image("./figures/spatial_patterns_normalized_flipped_eff_tirs.tif")
plot(im,axes=FALSE)
mtext(side = 3, line = -1.2, adj = 0.04, 'A', font = 2)
im <- load.image("./figures/spatial_patterns_normalized_flipped_eff_ulv.tif")
plot(im,axes=FALSE)
mtext(side = 3, line = -1.2, adj = 0.04, 'B', font = 2)
##par(mar = c(0,0.1,1.3,0.1))
dev.off()
