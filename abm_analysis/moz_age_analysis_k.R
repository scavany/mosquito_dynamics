library(moments)
load("age_dists.RData")
moments.arr <- list()

for (name in names(age.sorted.list)) {
    print(name)
    moments.arr[[name]] <- array(dim=c(nrow(age.sorted.list[[1]]),4))
    for (ii in 1:nrow(age.sorted.list[[name]])) {
        print(ii)
        temp.age.vec <- rep(as.numeric(names(age.sorted.list[[name]][ii,])),
                            as.numeric(age.sorted.list[[name]][ii,]))
        moments.arr[[name]][ii,1] <- mean(temp.age.vec)
        moments.arr[[name]][ii,2] <- sd(temp.age.vec)
        moments.arr[[name]][ii,3] <- skewness(temp.age.vec)
        moments.arr[[name]][ii,4] <- kurtosis(temp.age.vec)
    }
}

start=100
## Also plot the bimodality coefficient https://en.wikipedia.org/wiki/Multimodal_distribution#Bimodality_coefficient
for (name in names(age.sorted.list)) {
    tiff(paste0("figures/normalized_age_distribution_moments_",
                name,".tif"),
         res=600, width=4152, height=4152, compression="lzw")
    par(mfrow=c(2,2))
    plot(2000+(start:nrow(age.sorted.list[[name]])-1)/365,
         moments.arr[[name]][start:4018,1], type='l',
         xlab="Year",ylab="Mean",bty="n")
    if(name!="none")abline(v=2002,col="red",lty="dashed")
    plot(2000+(start:nrow(age.sorted.list[[name]])-1)/365,
         moments.arr[[name]][start:4018,2], type='l',
         xlab="Year",ylab="Standard deviation",bty="n")
    if(name!="none")abline(v=2002,col="red",lty="dashed")
    plot(2000+(start:nrow(age.sorted.list[[name]])-1)/365,
         moments.arr[[name]][start:4018,3], type='l',
         xlab="Year",ylab="Skewness",bty="n")
    if(name!="none")abline(v=2002,col="red",lty="dashed")
    plot(2000+(start:nrow(age.sorted.list[[name]])-1)/365,
         moments.arr[[name]][start:4018,4], type='l',
         xlab="Year",ylab="Kurtosis",bty="n")
    if(name!="none")abline(v=2002,col="red",lty="dashed")
    dev.off()

    tiff(paste0("figures/normalized_age_distribution_bimodality_",
                name,".tif"),
         res=600, width=4152, height=4152, compression="lzw")
    plot(2000+(start:nrow(age.sorted.list[[name]])-1)/365,
         (moments.arr[[name]][start:4018,3]^2 + 1)/moments.arr[[name]][start:4018,4], type='l',
         xlab="Year",ylab="Bimodality coefficient",bty="n")
    if(name!="none")abline(v=2002,col="red",lty="dashed")
    dev.off()
}
