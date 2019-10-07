## IRS, ULV, None age-dist bar charts
ages <- read.csv("moz_age_dist_het/sim_00000_moz_ages.csv")
age.sorted <- t(ages[order(ages$Age),])
colnames(age.sorted) <- age.sorted["Age",]
age.sorted <- age.sorted[-1,]

days.display <- seq(1, nrow(age.sorted), 1)
i <- 0
for (day in rownames(age.sorted)[days.display]) {
    age.sorted.plot <- age.sorted[day,]/sum(age.sorted[day,])
    jpeg(paste0("age_dist_", i, ".jpg"), quality=100)
    barplot(age.sorted.plot[1:60], ylim=c(0,0.3))
    if (i>=730 & i<770) {
        text(40, 0.2, "SPRAYING", col="red")
    }
    dev.off()
    i <- i+1
}

age.sorted.sum <- apply(age.sorted, 1, sum)
age.sorted.norm <- age.sorted/age.sorted.sum
image(age.sorted[100:4000,1:30])







##### needs expanding to do each one, and differences. Zone-level heterogeneity.
setwd("~/Desktop")
spatial.patterns <- array(dim=c(41, 35)) 
for (i in seq(41)) {
    day.num <- toString((i-1)*100)
    temp.data <- read.csv(paste0("moz_age_dist_het/heterogeneity_sim_00020_day_", day.num, ".csv"))
    temp.data$zoneID <- as.numeric(substring(temp.data$locID, 1, 2))
    for (j in seq(35)) {
        spatial.patterns[i,j] <- sum(temp.data$mozzes[which(temp.data$zoneID==j)])
    }
}

dates <- seq(as.Date("2000-01-01"), as.Date("2010-12-31"), "1 day")[seq(1,4001,100)]

image(x=dates, y=seq(1,35), z=spatial.patterns, col=plasma(200))

spatial.patterns.sum <- apply(spatial.patterns, 1, sum)
spatial.patterns.normalized <- spatial.patterns/spatial.patterns.sum
image(x=dates, y=seq(1,35), z=log(spatial.patterns.normalized), col=plasma(200), ylab="Zone", xlab="Year")


##k-functions
library(spatstat)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
setwd("~/Desktop/moz_age_dist_het")
xyz <- read.csv("heterogeneity_sim_00020_day_1000.csv")
iquitos <- readOGR("/home/sean/Documents/zika_project/shape_files_2/zonas35_ge.shp")
iquitos.ppp <- ppp(xyz$xCor, xyz$yCor, window=as.owin(iquitos))
marks(iquitos.ppp) <- xyz$mozzes
plot(Kest(iquitos.ppp))
plot(Kmulti(iquitos.ppp, rep(T, npoints(iquitos.ppp)),rep(T, npoints(iquitos.ppp))))


##### Clustering Gs staistic
library(akima)
library(rgdal)
library(spdep)
library(sp)
library(viridis)
library(mgcv)
setwd("~/Desktop/moz_age_dist_het")
xyz <- read.csv("heterogeneity_sim_00000_day_900.csv")
xycoords <- cbind(xyz$xCor, xyz$yCor)
nb30 <- dnearneigh(xycoords, 0, 30)
nb120 <- dnearneigh(xycoords, 0, 120)
G120 <- localG(xyz$mozzes, nb2listw(nb120, style="B"), )
G120s <- localG(log(xyz$mozzes+1), nb2listw(include.self(nb120),
                                  style="B"))
G30s <- localG(xyz$mozzes, nb2listw(include.self(nb30),
                                  style="B"))

iquitos <- readOGR("/home/sean/Documents/zika_project/shape_files_2/zonas35_ge.shp")
temp.grid <- makegrid(iquitos, n=100000)
iquitos.grid.full <- SpatialPoints(temp.grid, proj4string = CRS(proj4string(iquitos)))
iquitos.grid <- iquitos.grid.full[iquitos,]
G120.interp <- interp(xyz$xCor, xyz$yCor, G120,
                      xo=unique(iquitos.grid.full$x1),
                      yo=unique(iquitos.grid.full$x2),
                      duplicate="strip")
G120s.interp <- interp(xyz$xCor, xyz$yCor, G120s,
                       xo=unique(iquitos.grid.full$x1),
                       yo=unique(iquitos.grid.full$x2),
                       duplicate="strip")
G30s.interp <- interp(xyz$xCor, xyz$yCor, G30s,
                      xo=unique(iquitos.grid.full$x1),
                      yo=unique(iquitos.grid.full$x2),
                      duplicate="strip")
moz.interp <- interp(xyz$xCor, xyz$yCor, xyz$mozzes,
                     xo=unique(iquitos.grid.full$x1),
                     yo=unique(iquitos.grid.full$x2),
                     duplicate="strip")
moz.gam <- gam(mozzes~s(xCor)+s(yCor), data=xyz)
moz.gam.grid <- predict(moz.gam,
                        newdata=list(xCor=rep(unique(iquitos.grid.full$x1),
                                              times=length(unique(iquitos.grid.full$x2))),
                                     yCor=rep(unique(iquitos.grid.full$x2),
                                              each=length(unique(iquitos.grid.full$x1)))))
moz.gam.grid.mat <- matrix(moz.gam.grid,
                           nrow=length(unique(iquitos.grid.full$x1)),
                           ncol=length(unique(iquitos.grid.full$x2)))
G120.bounded <- G120.interp
moz.bounded <- moz.interp
G120s.bounded <- G120s.interp
G30s.bounded <- G30s.interp
for (i in seq(length(moz.interp$x))) {
    print(i)
    x.locs <- which(iquitos.grid$x1==moz.interp$x[i])
    y.vals <- iquitos.grid$x2[x.locs]
    for (j in seq(length(moz.interp$y))) {      
        if (!(moz.interp$y[j] %in% y.vals)) {
            print(paste("making a change", i, j))
            moz.gam.grid.mat[i,j] <- NA
            G120.bounded$z[i,j] <- NA
            G120s.bounded$z[i,j] <- NA
            G30s.bounded$z[i,j] <- NA
            moz.bounded$z[i,j] <- NA
        }
    }
}

image(G120.bounded, col=plasma(100))
plot(iquitos, add=T)
tiff(filename = "~/Desktop/moz_age_dist_het/G120s.tif",
          width = 4000, height = 4000,
          compression = "lzw", res = 600)
image(G120s.bounded, col=plasma(100), axes=F)
plot(iquitos, add=T)
dev.off()
image(G30s.bounded, col=plasma(100))
plot(iquitos, add=T)
image(moz.bounded, col=plasma(100))
plot(iquitos, add=T)
image(unique(iquitos.grid.full$x1),unique(iquitos.grid.full$x2),
      moz.gam.grid.mat, col=plasma(100))
plot(iquitos, add=T)
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
    image(G120s.bounded, col=plasma(100), axes=F)
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
image(G120s.bounded, col=plasma(100), axes=F)
plot(iquitos, add=T)
dev.off()

