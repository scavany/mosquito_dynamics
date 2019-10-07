setwd("~/Desktop")
imm.files <- list.files(pattern="immature.csv")
spr.files <- list.files(pattern="spray.csv")
n.files <- 1000

parms <- read.csv("parameters.csv")

n.timepts <- dim(read.csv(imm.files[1]))[1]
output <- array(dim=c(n.timepts, 2, n.files))

for (i in seq(length(imm.files))) {
    loc <- as.numeric(substring(imm.files[i], 5,9))
    output[,1,loc+1] <- read.csv(imm.files[i])[,"Adults"]
    output[,2,loc+1] <- read.csv(spr.files[i])[,2]
}

gam.moz <- read.csv("../Documents/zika_project/mosquito_dynamics/output/mosquito_ts.csv")[,1]
tvec <- seq(n.timepts)

plot(tvec, gam.moz[1:n.timepts],
     ylim=c(0, max(c(max(output[,1,], na.rm=T), gam.moz))),
     type='l', col='red')

for (i in seq(length(imm.files))) {
    loc <- as.numeric(substring(imm.files[i], 5,9))
    lines(tvec, output[,1,i], col="grey")
}

first.spray <- apply(output[,2,], 2, function(x) which(x>0)[1])
second.spray <- apply(output[366:n.timepts,2,], 2, function(x) 365+which(x>0)[1])

pre.abundance <- output[2,1,]
post.abundance <- apply(output[3:63, 1, ], 2, min)
drop.abundance <- (pre.abundance-post.abundance)/pre.abundance
nearest <- which.min(abs(drop.abundance-0.6))

parms[nearest, ]

pre.abundance.2 <- output[366,1,]
post.abundance.2 <- apply(output[367:437, 1, ], 2, min)
drop.abundance.2 <- 1 - post.abundance.2/pre.abundance.2
nearest.2 <- which.min(abs(drop.abundance.2-0.6))

parms[nearest.2, ]

drop.abundance.2[nearest]
drop.abundance[nearest.2]

nearest.all <- which.min(abs(drop.abundance.2-0.6)+abs(drop.abundance-0.6))
