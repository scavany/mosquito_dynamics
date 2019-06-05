#=================================================================
# User's input
#=================================================================
rm(list = ls())
library(here)
setwd(here())
setwd("/home/sean/Documents/zika_project/mosquito_dynamics/scripts/")

source('response_curves.R')

fieldcorxn <- 0.08969818

#=================================================================
# Compute the entomological parameters based on temperature 
#=================================================================
# Load emergence time-series: Emergence_TS from generate_emergence_file.R
# The rest of the variables should adjust to the emergence time-series
offsets = c("1.1","1.35","1.6","1.7","1.9", "2")
o_f = offsets[2]
#eval(parse(text = sprintf("load(\"../data/Emergence%s%s\")",o_f,".RData")))
# load('../data/Emergence.RData')
load('time_dependent_series.RData')

Iquitos.climate = read.csv('../data/Iquitos_Climate_Bobby.csv')
# startIndex = which(Iquitos.climate$date == "1999-02-28") # Assuming it starts in that date, change if necessary
startIndex = which(Iquitos.climate$date == "2000-01-01") 
temperature = Iquitos.climate$temperature_mean[startIndex:(startIndex + length(larval_development_rate)-1) ]
temp.mean <- mean(temperature,na.rm = TRUE)
temperature[is.na(temperature)] <- temp.mean

eip.mean = array(data = 0, dim = c(length(larval_development_rate),1))
eip.mu = array(data = 0, dim = c(length(larval_development_rate),1))
bite.first = array(data = 0, dim = c(length(larval_development_rate),1))
bite.second = array(data = 0, dim = c(length(larval_development_rate),1))
moz.death = array(data = 0, dim = c(length(larval_development_rate),1))

for (t in 1:length(larval_development_rate)){
  eip.mu[t] = log(eip(temperature[t], 402, .1)$mu)
  eip.mean[t] = eip(temperature[t], 402, .1)$mean
  bite.first[t] = biterate.1st(temperature[t])
  bite.second[t] = biterate.2nd(temperature[t])
}

EIP.data = data.frame(mu = eip.mu)
mosquito.data = data.frame(EIP = eip.mu, firstBite = bite.first, secondBite = bite.second, eggDeath = egg_mortality, larvaeDeath = larval_mortality, pupalDeath = pupal_mortality,  death = DeathRate, eggDev=egg_development_rate, larvalDev = larval_development_rate, pupalDev = pupal_development_rate, adultDev=gonotrophic_cycle_rate, larvalCap=larval_capacity_timeseries)
#mosquito.data = data.frame(EIP = eip.mu, firstBite = bite.first, secondBite = bite.second, eggDeath = egg_mortality, larvaeDeath = larval_mortality, pupalDeath = pupal_mortality,  death = DeathRate, eggDev=egg_development_timeseries_nn6, larvalDev = larval_development_rate, pupalDev = pupal_development_rate, adultDev=gonotrophic_cycle_rate)
mosquito.data <- rbind(mosquito.data, colMeans(mosquito.data))
write.csv(mosquito.data, sprintf("../output/AegyptiPars_%s_offset.csv",o_f),row.names = FALSE, quote = FALSE)

