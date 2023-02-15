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
sim.nums <- c(tirs="00000", ulv="00010", none="00020")
sim.names <- names(sim.nums)

## IRS, ULV, None age-dist bar charts
age.sorted.list <- list()
cols=viridis(4)
cols.alpha=viridis(4, alpha=0.33)
for (sim.num in sim.nums) { 
    ages <- read.csv(paste0("moz_age_dist_het/sim_",sim.num,"_moz_ages.csv"))
    mozzes <- read.csv(paste0("moz_age_dist_het/sim_",sim.num,"_immature.csv"))[,"Adults"]
    ylims <- c(0,7e5)
    age.sorted <- t(ages[order(ages$Age),])
    colnames(age.sorted) <- age.sorted["Age",]
    age.sorted <- age.sorted[-1,]
    age.sorted.list[[names(which(sim.nums==sim.num))]] <- age.sorted
    days.display <- seq(1, nrow(age.sorted), 1)

    for (day in rownames(age.sorted)[days.display]) {
        i <- as.numeric(substr(day, start=5, stop=nchar(day)))
        age.sorted.plot <- age.sorted[day,]/sum(age.sorted[day,])
        png(paste0("figures/sim_", sim.num, "_age_dist_", i, ".png"),
            height=600, width=600)
        par(fig=c(0,1,0,1))
        barplot(age.sorted.plot[1:60], ylim=c(0,0.3), col=cols[1],
                xlab="Age (days)", ylab="Proportion")
        if (names(which(sim.nums==sim.num)) == "tirs") {
            if (i>=731 & i<=770) {
                text(15, 0.25, "SPRAYING",
                     col=cols[2])
            }
        } else if (names(which(sim.nums==sim.num)) == "ulv") {
            if (i>=731 & i<=758) {
                text(15, 0.2, "SPRAYING",
                     col=cols[2])
            }
        }
        par(fig = c(0.5,0.95, 0.5, 0.95), new = T)  
        plot(2000+days.display/365, mozzes[days.display],
             xlab="Year", ylab="Abundance",
             type='l', col=cols[3])
        if (names(which(sim.nums==sim.num)) == "tirs") {
            timepts <- 2000+seq(731, 770)/365
            ## polygon(c(timepts, rev(timepts)),
            ##         c(rep(ylims[1],length(timepts)), rep(ylims[2], length(timepts))),
            ##         col=cols.alpha[2], border=F)
        } else if (names(which(sim.nums==sim.num)) == "ulv") {
            timepts <- 2000+seq(731, 770)/365
            ## polygon(c(timepts, rev(timepts)),
            ##         c(rep(ylims[1],length(timepts)), rep(ylims[2], length(timepts))),
            ##         col=cols.alpha[2], border=F)
        }
        abline(v=2000+i/365, lty="dashed")
        dev.off()
    }
    system(paste0("ffmpeg -i figures/sim_",
                  sim.num,
                  "_age_dist_%d.png -c:v libx264 -y figures/sim_",
                  sim.num, "_age_dist.mp4"))
}

age.sorted.sum <- lapply(age.sorted.list, function(x) apply(x, 1, sum))
age.sorted.norm <- Map("/",age.sorted.list, age.sorted.sum)

save(age.sorted.list, age.sorted.norm, file="age_dists.RData")

load("age_dists.RData")
##image.plot(2000+100:4000/365, 1:30-1, age.sorted[100:4000,1:30], col=plasma(100),
##      xlab="Year", ylab="Age (days)")
for (i in seq(length(age.sorted.list))) {
    tiff(paste0("figures/normalized_age_distribution_",
                names(age.sorted.list)[[i]],".tif"),
         res=600, width=4152, height=4152, compression="lzw")
    image.plot(2000+100:4000/365, 1:30-1, age.sorted.norm[[i]][100:4000,1:30], col=plasma(100),
               xlab="Year", ylab="Age (days)",
               main=paste0("Proportion of female adults at each age\n",
                           ifelse(names(age.sorted.list)[[i]] != "none",paste("following",toupper(names(age.sorted.list)[[i]])),"")))
    if (names(age.sorted.list)[[i]] != "none") {
        abline(v=2002,col="red",lwd=2,lty="dashed")
    }
    dev.off()
    if (names(age.sorted.list)[[i]] != "none") {
        tiff(paste0("figures/normalized_age_difference_",
                    names(age.sorted.list)[[i]],".tif"),
             res=600, width=4152, height=4152, compression="lzw")
        image.plot(2000+100:4000/365, 1:30-1,
              age.sorted.norm[[i]][100:4000,1:30] - age.sorted.norm[["none"]][100:4000,1:30],
              col=plasma(100),
              xlab="Year", ylab="Age (days)",
              main=paste0("Difference in proportion of female adults in each zone\n",
                          "following ",toupper(names(age.sorted.list)[[i]]),""))
        abline(v=2002,col="red",lwd=2,lty="dashed")
        dev.off()
    }
}


##### needs expanding to do each one, and differences. Zone-level heterogeneity.
setwd("~/Documents/zika_project/mosquito_dynamics/abm_analysis")
spatial.patterns <- array(dim=c(3, 41, 35)) 
dates <- seq(as.Date("2000-01-01"), as.Date("2010-12-31"), "1 day")[seq(1,4018,100)]
for (sim.num in sim.nums) {
    k <- which(sim.nums==sim.num)
    for (i in seq(41)) {
        day.num <- toString((i-1)*100)
        temp.data <- read.csv(paste0("moz_age_dist_het/heterogeneity_sim_",
                                     sim.num, "_day_", day.num, ".csv"))
        temp.data$zoneID <- as.numeric(substring(temp.data$locID, 1, 2))
        for (j in seq(35)) {
            spatial.patterns[k,i,j] <- sum(temp.data$mozzes[which(temp.data$zoneID==j)])
        }
    }
    
}

save(spatial.patterns, file="spatial_patterns.RData")

load("spatial_patterns.RData")
## ##map and pattern baseline
## sim.num <- "00020"
## k <- 3
## tiff(paste0("figures/spatial_patterns_baseline_summary.tif"),
##      res=600, width=4152, height=2076, compression="lzw")
## par(mfrow=c(1,2))

## spatial.patterns.sum <- apply(spatial.patterns[k,,], 1, sum)
## spatial.patterns.normalized <- spatial.patterns[k,,]/spatial.patterns.sum
## image.plot(x=dates, y=seq(1,35), z=spatial.patterns.normalized,
##       col=plasma(200), ylab="Zone", xlab="Year")

for (sim.num in sim.nums) {
    k <- which(sim.nums==sim.num)
    tiff(paste0("figures/spatial_patterns_",
                names(sim.nums)[k],".tif"),
         res=600, width=4152, height=4152, compression="lzw")
    image.plot(x=dates, y=seq(1,35), z=spatial.patterns[k,,],
               col=plasma(200), ylab="Zone", xlab="Year",
               main="Total female adult abundance")
    if (names(which(sim.nums==sim.num)) != "none") {
        abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
    }
    dev.off()
    spatial.patterns.sum <- apply(spatial.patterns[k,,], 1, sum)
    spatial.patterns.normalized <- spatial.patterns[k,,]/spatial.patterns.sum
    tiff(paste0("figures/spatial_patterns_normalized_",
                names(sim.nums)[k],".tif"),
         res=600, width=4152, height=4152, compression="lzw")
    image.plot(x=dates, y=seq(1,35), z=spatial.patterns.normalized,
          col=plasma(200), ylab="Zone", xlab="Year",zlim=c(0,max(spatial.patterns.normalized)),
          main=paste0("Proportion of female adults in each zone\n",
                      ifelse(names(which(sim.nums==sim.num)) != "none",paste("following",toupper(names(sim.nums)[k])),"")))
    if (names(which(sim.nums==sim.num)) != "none") {
        abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
    }
    dev.off()
    if (names(which(sim.nums==sim.num)) == "tirs") {
        tiff(paste0("figures/spatial_patterns_normalized_",
                    names(sim.nums)[k],"_ulv_colorbar.tif"),
             res=600, width=4152, height=4152, compression="lzw")
        k.ulv <- which(names(sim.nums)=="ulv")
        spatial.patterns.sum.ulv <- apply(spatial.patterns[k.ulv,,], 1, sum)
        spatial.patterns.normalized.ulv <- spatial.patterns[k.ulv,,]/spatial.patterns.sum.ulv
        image.plot(x=dates, y=seq(1,35), z=spatial.patterns.normalized,
                   col=plasma(200), ylab="Zone", xlab="Year",zlim=c(0,max(spatial.patterns.normalized.ulv)),
                   main=paste0("Proportion of female adults in each zone\n",
                               ifelse(names(which(sim.nums==sim.num)) != "none",paste("following",toupper(names(sim.nums)[k])),"")))
        abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
        dev.off()
    }
    
    tiff(paste0("figures/spatial_patterns_normalized_log_",
                names(sim.nums)[k],".tif"),
         res=600, width=4152, height=4152, compression="lzw")
    image(x=dates, y=seq(1,35), z=log(spatial.patterns.normalized),
          col=plasma(200), ylab="Zone", xlab="Year",
          main="Proportion of female adults in each zone, log-scale")
    if (names(which(sim.nums==sim.num)) != "none") {
        abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
    }
    dev.off()
    if (names(which(sim.nums==sim.num)) != "none") {
        k0 <- which(sim.names=="none")
        spatial.patterns.sum.0 <- apply(spatial.patterns[k0,,], 1, sum)
        spatial.patterns.normalized.0 <- spatial.patterns[k0,,]/spatial.patterns.sum.0
        tiff(paste0("figures/spatial_patterns_difference_",
                    names(sim.nums)[k],".tif"),
             res=600, width=4152, height=4152, compression="lzw")
        image.plot(x=dates, y=seq(1,35),
                   z=spatial.patterns[k0,,]-spatial.patterns[k,,],
                   col=plasma(200), ylab="Zone", xlab="Year",
                   main="Change in abundance of female adults following spraying")
        if (names(which(sim.nums==sim.num)) != "none") {
            abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
        }
        dev.off()
        tiff(paste0("figures/spatial_patterns_normalized_difference_",
                    names(sim.nums)[k],".tif"),
             res=600, width=4152, height=4152, compression="lzw")
        image.plot(x=dates, y=seq(1,35),
              z=spatial.patterns.normalized.0-spatial.patterns.normalized,
              col=plasma(200), ylab="Zone", xlab="Year",
              main="Change in proportion of female adults following spraying")
        if (names(which(sim.nums==sim.num)) != "none") {
            abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
        }
        dev.off()
        tiff(paste0("figures/spatial_patterns_normalized_difference_log_",
                    names(sim.nums)[k],".tif"),
             res=600, width=4152, height=4152, compression="lzw")
        image(x=dates, y=seq(1,35),
              z=log(abs(spatial.patterns.normalized.0-spatial.patterns.normalized)),
              col=plasma(200), ylab="Zone", xlab="Year",
              main="Change in proportion of female adults following spraying,\nlog-scale")
        if (names(which(sim.nums==sim.num)) != "none") {
            abline(v=as.Date("2002-01-01"),col="red",lwd=2,lty="dashed")
        }
        dev.off()
    }
}

## ##k-functions
## setwd("~/Desktop/moz_age_dist_het")
## xyz <- read.csv("heterogeneity_sim_00020_day_1000.csv")
## iquitos <- readOGR("/home/sean/Documents/zika_project/shape_files_2/zonas35_ge.shp")
## iquitos.ppp <- ppp(xyz$xCor, xyz$yCor, window=as.owin(iquitos))
## marks(iquitos.ppp) <- xyz$mozzes
## plot(Kest(iquitos.ppp))
## plot(Kmulti(iquitos.ppp, rep(T, npoints(iquitos.ppp)),rep(T, npoints(iquitos.ppp))))


##### Mapping abundance and (deprecated) clustering statistic
setwd("~/Documents/zika_project/mosquito_dynamics/abm_analysis/moz_age_dist_het")
map.days <- seq(0,4018,100)
dates <- seq(as.Date("2000-01-01"), as.Date("2010-12-31"), "1 day")[seq(1,4018,100)]
xyz <- read.csv("heterogeneity_sim_00000_day_0.csv")
xycoords <- cbind(xyz$xCor, xyz$yCor)
z.lims <- c(0,17)
## nb30 <- dnearneigh(xycoords, 0, 30)
## nb120 <- dnearneigh(xycoords, 0, 120)
## G120 <- localG(xyz$mozzes, nb2listw(nb120, style="B"), )
## G120s <- localG(log(xyz$mozzes+1), nb2listw(include.self(nb120),
##                                   style="B"))
## G30s <- localG(xyz$mozzes, nb2listw(include.self(nb30),
##                                   style="B"))

iquitos <- readOGR("/home/sean/Documents/zika_project/shape_files_2/zonas35_ge.shp")
temp.grid <- makegrid(iquitos, n=100000)
iquitos.grid.full <- SpatialPoints(temp.grid, proj4string = CRS(proj4string(iquitos)))
iquitos.grid <- iquitos.grid.full[iquitos,]
## G120.interp <- interp(xyz$xCor, xyz$yCor, G120,
##                       xo=unique(iquitos.grid.full$x1),
##                       yo=unique(iquitos.grid.full$x2),
##                       duplicate="strip")
## G120s.interp <- interp(xyz$xCor, xyz$yCor, G120s,
##                        xo=unique(iquitos.grid.full$x1),
##                        yo=unique(iquitos.grid.full$x2),
##                        duplicate="strip")
## G30s.interp <- interp(xyz$xCor, xyz$yCor, G30s,
##                       xo=unique(iquitos.grid.full$x1),
##                       yo=unique(iquitos.grid.full$x2),
##                       duplicate="strip")
moz.interp <- interp(xyz$xCor, xyz$yCor, xyz$mozzes,
                     xo=unique(iquitos.grid.full$x1),
                     yo=unique(iquitos.grid.full$x2),
                     duplicate="strip")
for (sim in names(sim.nums)) {
    print(sim)
    for (jjj in seq(length(map.days))) {
        print(jjj)
    ## ## if (map.days$When[jjj] == "before") {
    ##     xyz <- read.csv(paste0("heterogeneity_sim_",sim.nums["none"],"_day_",
    ##                            map.days[jjj], ".csv"))
    ##     moz.gam <- gam(mozzes~s(xCor)+s(yCor), data=xyz)
    ##     moz.gam.grid <- predict(moz.gam,
    ##                             newdata=list(xCor=rep(unique(iquitos.grid.full$x1),
    ##                                                   times=length(unique(iquitos.grid.full$x2))),
    ##                                          yCor=rep(unique(iquitos.grid.full$x2),
    ##                                                   each=length(unique(iquitos.grid.full$x1)))))
    ##     moz.gam.grid.mat <- matrix(moz.gam.grid,
    ##                                nrow=length(unique(iquitos.grid.full$x1)),
    ##                                ncol=length(unique(iquitos.grid.full$x2)))
    ##     ## G120.bounded <- G120.interp
    ##     ## moz.bounded <- moz.interp
    ##     ## G120s.bounded <- G120s.interp
    ##     ## G30s.bounded <- G30s.interp
    ##     for (i in seq(length(moz.interp$x))) {
    ##         print(i)
    ##         x.locs <- which(iquitos.grid$x1==moz.interp$x[i])
    ##         y.vals <- iquitos.grid$x2[x.locs]
    ##         for (j in seq(length(moz.interp$y))) {      
    ##             if (!(moz.interp$y[j] %in% y.vals)) {
    ##                 ## print(paste("making a change", i, j))
    ##                 moz.gam.grid.mat[i,j] <- NA
    ##                 ## G120.bounded$z[i,j] <- NA
    ##                 ## G120s.bounded$z[i,j] <- NA
    ##                 ## G30s.bounded$z[i,j] <- NA
    ##                 ## moz.bounded$z[i,j] <- NA
    ##             }
    ##         }
    ##     }

    ##     ## image.plot(G120.bounded, col=plasma(100))
    ##     ## plot(iquitos, add=T)
    ##     ## tiff(filename = "~/Desktop/moz_age_dist_het/G120s.tif",
    ##     ##           width = 4000, height = 4000,
    ##     ##           compression = "lzw", res = 600)
    ##     ## image.plot(G120s.bounded, col=plasma(100), axes=F)
    ##     ## plot(iquitos, add=T)
    ##     ## dev.off()
    ##     ## image.plot(G30s.bounded, col=plasma(100))
    ##     ## plot(iquitos, add=T)
    ##     ## image.plot(moz.bounded, col=plasma(100))
    ##     ## plot(iquitos, add=T)
    ##     png(paste0("../figures/moz_density_", map.days$When[jjj], ".png"),
    ##         height=600, width=600)
    ##     image.plot(unique(iquitos.grid.full$x1),unique(iquitos.grid.full$x2),
    ##           moz.gam.grid.mat, col=plasma(100), axes=F, xlab=NA, ylab=NA)
    ##     plot(iquitos, add=T)
    ##     dev.off()
    ## ## } else {
        xyz <- read.csv(paste0("heterogeneity_sim_",sim.nums[sim],"_day_",
                               map.days[jjj], ".csv"))
        moz.gam <- gam(mozzes~s(xCor)+s(yCor), data=xyz)
        moz.gam.grid <- predict(moz.gam,
                                newdata=list(xCor=rep(unique(iquitos.grid.full$x1),
                                                      times=length(unique(iquitos.grid.full$x2))),
                                             yCor=rep(unique(iquitos.grid.full$x2),
                                                      each=length(unique(iquitos.grid.full$x1)))))
        moz.gam.grid.mat <- matrix(moz.gam.grid,
                                   nrow=length(unique(iquitos.grid.full$x1)),
                                   ncol=length(unique(iquitos.grid.full$x2)))
        ## G120.bounded <- G120.interp
        ## moz.bounded <- moz.interp
        ## G120s.bounded <- G120s.interp
        ## G30s.bounded <- G30s.interp
        for (i in seq(length(moz.interp$x))) {
            ## print(i)
            x.locs <- which(iquitos.grid$x1==moz.interp$x[i])
            y.vals <- iquitos.grid$x2[x.locs]
            for (j in seq(length(moz.interp$y))) {      
                if (!(moz.interp$y[j] %in% y.vals)) {
                    ## print(paste("making a change", i, j))
                    moz.gam.grid.mat[i,j] <- NA
                    ## G120.bounded$z[i,j] <- NA
                    ## G120s.bounded$z[i,j] <- NA
                    ## G30s.bounded$z[i,j] <- NA
                    ## moz.bounded$z[i,j] <- NA
                }
            }
        }
        
        ## image.plot(G120.bounded, col=plasma(100))
        ## plot(iquitos, add=T)
        ## tiff(filename = "~/Desktop/moz_age_dist_het/G120s.tif",
        ##           width = 4000, height = 4000,
        ##           compression = "lzw", res = 600)
        ## image.plot(G120s.bounded, col=plasma(100), axes=F)
        ## plot(iquitos, add=T)
        ## dev.off()
        ## image.plot(G30s.bounded, col=plasma(100))
        ## plot(iquitos, add=T)
        ## image.plot(moz.bounded, col=plasma(100))
        ## plot(iquitos, add=T)
        png(paste0("../figures/moz_density_", jjj,
                   "_", sim, ".png"),
            height=600, width=600)
        image.plot(unique(iquitos.grid.full$x1),unique(iquitos.grid.full$x2),
              moz.gam.grid.mat, col=plasma(10*diff(z.lims)), zlim=z.lims,
              axes=F, xlab=NA, ylab=NA)
        plot(iquitos, add=T)
        text(min(unique(iquitos.grid.full$x1)),max(unique(iquitos.grid.full$x2)),
             dates[jjj], adj=c(0,1))
        dev.off()
    }
    ## }
    system(paste0("ffmpeg -r 1 -i ../figures/moz_density_%d_",sim,
                  ".png -c:v libx264 -y ../figures/sim_",
                  sim, "_moz_density.mp4"))
}

##### Mapping normalized abundance
setwd("~/Documents/zika_project/mosquito_dynamics/abm_analysis/moz_age_dist_het")
map.days <- seq(0,4018,100)
dates <- seq(as.Date("2000-01-01"), as.Date("2010-12-31"), "1 day")[seq(1,4018,100)]
xyz <- read.csv("heterogeneity_sim_00000_day_0.csv")
xycoords <- cbind(xyz$xCor, xyz$yCor)
z.lims <- c(0,17)

iquitos <- readOGR("/home/sean/Documents/zika_project/shape_files_2/zonas35_ge.shp")
temp.grid <- makegrid(iquitos, n=100000)
iquitos.grid.full <- SpatialPoints(temp.grid, proj4string = CRS(proj4string(iquitos)))
iquitos.grid <- iquitos.grid.full[iquitos,]
moz.interp <- interp(xyz$xCor, xyz$yCor, xyz$mozzes,
                     xo=unique(iquitos.grid.full$x1),
                     yo=unique(iquitos.grid.full$x2),
                     duplicate="strip")
for (sim in names(sim.nums)) {
    print(sim)
    for (jjj in seq(length(map.days))) {
        print(jjj)
        xyz <- read.csv(paste0("heterogeneity_sim_",sim.nums[sim],"_day_",
                               map.days[jjj], ".csv"))
        moz.gam <- gam(mozzes~s(xCor)+s(yCor), data=xyz)
        moz.gam.grid <- predict(moz.gam,
                                newdata=list(xCor=rep(unique(iquitos.grid.full$x1),
                                                      times=length(unique(iquitos.grid.full$x2))),
                                             yCor=rep(unique(iquitos.grid.full$x2),
                                                      each=length(unique(iquitos.grid.full$x1)))))
        moz.gam.grid.mat <- matrix(moz.gam.grid,
                                   nrow=length(unique(iquitos.grid.full$x1)),
                                   ncol=length(unique(iquitos.grid.full$x2)))
        for (i in seq(length(moz.interp$x))) {
            x.locs <- which(iquitos.grid$x1==moz.interp$x[i])
            y.vals <- iquitos.grid$x2[x.locs]
            for (j in seq(length(moz.interp$y))) {      
                if (!(moz.interp$y[j] %in% y.vals)) {
                    moz.gam.grid.mat[i,j] <- NA
                }
            }
        }
        png(paste0("../figures/moz_density_", jjj,
                   "_", sim, "_normalized.png"),
            height=600, width=600)
        image.plot(unique(iquitos.grid.full$x1),unique(iquitos.grid.full$x2),
                   moz.gam.grid.mat, col=plasma(10*diff(z.lims)),
                   axes=F, xlab=NA, ylab=NA,
                   main=paste0("Female adult abundance on\n",dates[jjj]))
        plot(iquitos, add=T)
        text(coordinates(iquitos), labels=seq(35), col="white")
        ## text(min(unique(iquitos.grid.full$x1)),max(unique(iquitos.grid.full$x2)),
        ##      dates[jjj], adj=c(0,1))
        dev.off()
    }
    system(paste0("ffmpeg -r 1 -i ../figures/moz_density_%d_",sim,
                  "_normalized.png -c:v libx264 -y ../figures/sim_",
                  sim, "_moz_density_normalized.mp4"))
}

## Do each timept for video
iquitos <- readOGR("/home/sean/Documents/zika_project/shape_files_2/zonas35_ge.shp")
temp.grid <- makegrid(iquitos, n=100000)
iquitos.grid.full <- SpatialPoints(temp.grid, proj4string = CRS(proj4string(iquitos)))
iquitos.grid <- iquitos.grid.full[iquitos,]
for (day in seq(0,4018,100)) {
    print(day)
    xyz <- read.csv(paste0("heterogeneity_sim_00020_day_",day,".csv"))
    xycoords <- cbind(xyz$xCor, xyz$yCor)
    nb120 <- dnearneigh(xycoords, 0, 120)
    G120s <- localG(log(xyz$mozzes+1), nb2listw(include.self(nb120),
                                                style="B"))
    G120s.interp <- interp(xyz$xCor, xyz$yCor, G120s,
                           xo=unique(iquitos.grid.full$x1),
                           yo=unique(iquitos.grid.full$x2),
                           duplicate="strip")
    G120s.bounded <- G120s.interp
    for (i in seq(length(G120.interp$x))) {
        x.locs <- which(iquitos.grid$x1==G120s.interp$x[i])
        y.vals <- iquitos.grid$x2[x.locs]
        for (j in seq(length(G120s.interp$y))) {      
            if (!(G120.interp$y[j] %in% y.vals)) {
                G120s.bounded$z[i,j] <- NA
            }
        }
    }

    tiff(filename = paste0("~/Desktop/moz_age_dist_het/G120s_",
                           day,".tif"),
          width = 4000, height = 4000,
         compression = "lzw", res = 600)
    image.plot(G120s.bounded, col=plasma(100), axes=F)
    plot(iquitos, add=T)
    dev.off()
}

##try IRS with day 800 and/or 900 - the ones following spraying.
setwd("~/Desktop/moz_age_dist_het/")
iquitos <- readOGR("/home/sean/Documents/zika_project/shape_files_2/zonas35_ge.shp")
temp.grid <- makegrid(iquitos, n=100000)
iquitos.grid.full <- SpatialPoints(temp.grid, proj4string = CRS(proj4string(iquitos)))
iquitos.grid <- iquitos.grid.full[iquitos,]
day <- 800
xyz <- read.csv(paste0("heterogeneity_sim_00000_day_",day,".csv"))
xycoords <- cbind(xyz$xCor, xyz$yCor)
nb120 <- dnearneigh(xycoords, 0, 120)
G120s <- localG(log(xyz$mozzes+1), nb2listw(include.self(nb120),
                                            style="B"))
G120s.interp <- interp(xyz$xCor, xyz$yCor, G120s,
                       xo=unique(iquitos.grid.full$x1),
                       yo=unique(iquitos.grid.full$x2),
                       duplicate="strip")
G120s.bounded <- G120s.interp
for (i in seq(length(G120.interp$x))) {
    x.locs <- which(iquitos.grid$x1==G120s.interp$x[i])
    y.vals <- iquitos.grid$x2[x.locs]
    for (j in seq(length(G120s.interp$y))) {      
        if (!(G120.interp$y[j] %in% y.vals)) {
            G120s.bounded$z[i,j] <- NA
            }
    }
}

tiff(filename = paste0("~/Desktop/moz_age_dist_het/IRS_G120s_",
                       day,".tif"),
     width = 4000, height = 4000,
     compression = "lzw", res = 600)
image.plot(G120s.bounded, col=plasma(100), axes=F)
plot(iquitos, add=T)
dev.off()

