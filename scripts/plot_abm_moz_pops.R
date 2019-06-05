setwd("/home/sean/Documents/zika_project/mosquito_dynamics/scripts/")

load("surfaces.RData")

data_egg = read.csv("../../outbreak_response/IquitoSim_code/output/sim_zonal_eggdev_foi.csv")
data_capacity = read.csv("../../outbreak_response/IquitoSim_code/output/sim_zonal0_foi.csv")
data_delay = read.csv("../../outbreak_response/IquitoSim_code/output/sim_zonal_foi.csv")

moz_egg = data_egg$MozSusceptible_1 + data_egg$MozExposed_1 +data_egg$MozInfectious_1
moz_capacity = data_capacity$MozSusceptible_1 + data_capacity$MozExposed_1 +data_capacity$MozInfectious_1
moz_delay = data_delay$MozSusceptible_1 + data_delay$MozExposed_1 +data_delay$MozInfectious_1

moz_egg_n = moz_egg/moz_egg[1]
moz_capacity_n = moz_capacity/moz_capacity[1]
moz_delay_n = moz_delay/moz_delay[1]

plot(adult_timeseries_trend[1:1000], ylim = c(0.95*min(c(moz_egg_n, moz_delay_n, moz_capacity_n, adult_timeseries_trend)),1.05*max(c(moz_egg_n, moz_delay_n, moz_capacity_n, adult_timeseries_trend))))
legend(0, y=8, legend=c("egg_development", "larval capacity", "delay DEs"), lty=c(1,1,1), col = c("red", "blue", "green"))
lines(moz_egg_n, col='red')
lines(moz_capacity_n, col='blue')
lines(moz_delay_n, col='green')


dev.copy(png,'mozplot.png')
dev.off()

plot(adult_timeseries_trend[1:1000], ylim = c(0.95*min(c(moz_delay_n, moz_capacity_n, adult_timeseries_trend)),1.05*max(c(moz_delay_n, moz_capacity_n, adult_timeseries_trend))))
legend(0, y=2.9, legend=c("larval capacity", "delay DEs"), lty=c(1,1,1), col = c("red", "blue", "green"))
lines(moz_capacity_n, col='blue')
lines(moz_delay_n, col='green')

dev.copy(png,'mozplot2.png')
dev.off()
