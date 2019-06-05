#=================================================================
# User's input
#=================================================================
rm(list = ls())
library(here)
setwd(here())
setwd("/home/sean/Documents/zika_project/mosquito_dynamics/scripts/")

source('response_curves.R')

#fieldcorxn <- 0.08969818

#=================================================================
# Compute the entomological parameters based on temperature 
#=================================================================
# Load emergence time-series: Emergence_TS from generate_emergence_file.R
# The rest of the variables should adjust to the emergence time-series
offsets = c("1.1","1.35","1.6","1.7","1.9", "2")
o_f = offsets[2]
#eval(parse(text = sprintf("load(\"../data/Emergence%s%s\")",o_f,".RData")))
# load('../data/Emergence.RData')
load('time_dependent_series_delay.RData')
DeathRate# = DeathRate[(burn_period+1):(length(DeathRate)-burn_period)]

Iquitos.climate = read.csv('../data/Iquitos_Climate_Bobby.csv')
# startIndex = which(Iquitos.climate$date == "1999-02-28") # Assuming it starts in that date, change if necessary
startIndex = which(Iquitos.climate$date == "2000-01-01") 
temperature = Iquitos.climate$temperature_mean[startIndex:(startIndex + length(DeathRate)-1) ]
temp.mean <- mean(temperature,na.rm = TRUE)
temperature[is.na(temperature)] <- temp.mean

eip.mean = array(data = 0, dim = c(length(DeathRate),1))
eip.mu = array(data = 0, dim = c(length(DeathRate),1))
bite.first = array(data = 0, dim = c(length(DeathRate),1))
bite.second = array(data = 0, dim = c(length(DeathRate),1))
moz.death = array(data = 0, dim = c(length(DeathRate),1))

for (t in 1:length(DeathRate)){
  eip.mu[t] = log(eip(temperature[t], 402, .1)$mu)
  eip.mean[t] = eip(temperature[t], 402, .1)$mean
  bite.first[t] = biterate.1st(temperature[t])
  bite.second[t] = biterate.2nd(temperature[t])
}

EIP.data = data.frame(mu = eip.mu)
mosquito.data = data.frame(EIP = eip.mu, firstBite = bite.first, secondBite = bite.second, death = DeathRate, emerge = emergence_1*correction_factor)
mosquito.data <- rbind(mosquito.data, colMeans(mosquito.data))
write.csv(mosquito.data, sprintf("../output/AegyptiPars_%s_offset_delay.csv",o_f),row.names = FALSE, quote = FALSE)

