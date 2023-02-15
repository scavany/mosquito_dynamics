#=============================================================================
# Author: Sean Cavany
# 09/15/2019 
# Produce figures analyzing ABM moz output
#
#=============================================================================
library(akima)
library(rgdal)
library(spdep)
library(sp)
library(viridis)
library(mgcv)
library(spatstat)
library(rgeos)
library(maptools)
library(fields)

setwd("~/Documents/zika_project/mosquito_dynamics/abm_analysis")
sim.num <- "00000"

spatial.patterns <- array(dim=c(41, 35)) 
dates <- seq(as.Date("2000-01-01"), as.Date("2010-12-31"), "1 day")[seq(1,4018,100)]
for (i in seq(41)) {
    day.num <- toString((i-1)*100)
    temp.data <- read.csv(paste0("moz_age_dist_het_tirs_increased_mozpop/heterogeneity_sim_",
                                 sim.num, "_day_", day.num, ".csv"))
    temp.data$zoneID <- as.numeric(substring(temp.data$locID, 1, 2))
    for (j in seq(35)) {
        spatial.patterns[i,j] <- sum(temp.data$mozzes[which(temp.data$zoneID==j)])
    }
}

## save(spatial.patterns, file="spatial_patterns_increased_mozpop.RData")

load("spatial_patterns_increased_mozpop.RData")

spatial.patterns.sum <- apply(spatial.patterns, 1, sum)
spatial.patterns.normalized <- spatial.patterns/spatial.patterns.sum
tiff("figures/spatial_patterns_normalized_tirs_increased_mozpop.tif",
     res=600, width=4152, height=4152, compression="lzw")
image.plot(x=dates, y=seq(1,35), z=spatial.patterns.normalized,
           col=plasma(200), ylab="Zone", xlab="Year")
abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
dev.off()
