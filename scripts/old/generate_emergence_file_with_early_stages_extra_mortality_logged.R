#======================================================================================
# Author: Guido Espana
# project: PLOS COMP BIO paper 2016 Iquitos 
#          GSK trial simulations
# Year: 2016
# 
# Code to derive the location and time specific emergence of mosquitoes in Iquitos, 
# based on weekly adult sampling data. 
#
# Files needed:
# ../data/locations_20140807.csv
# ../data/Iquitos_Climate_Bobby.csv
# response_curves.R
# algam_85re.Rdata
# EstimateSurface[1.1,1.35,1.6,1.7,1.9.2].RData
# 
#======================================================================================
# Load libraries and source files --------
#======================================================================================
rm(list = ls())

# Set your working directory to the scripts folder
library(here)
setwd(here())

setwd("/home/sean/Documents/zika_project/mosquito_dynamics/scripts/")

require(mgcv)
require(scam)
require(deSolve)
require(pspline)
library(doParallel)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dlnm,
               mgcv,
               scam,
               mapproj,
               maptools,
               spatstat,
               RgoogleMaps,
               RColorBrewer,
               plotrix,
               rgdal,
               rgeos,
               ggplot2,
               raster,
               gridExtra,
               splines,
               zoo)

source('./generate_emergence_functions_with_early_stages.R')
source('response_curves.R')


#======================================================================================
# Load data -----------------
#======================================================================================
dvec = seq(as.Date("2000-01-01"),as.Date("2010-12-31"),by="1 day")
tvec = seq(as.Date("2000-01-01"),as.Date("2010-12-31"),by="1 day")
# tvec = seq(as.POSIXct("2000-01-01"),as.POSIXct("2010-12-31 23:00"),by="hours")
initial_date = as.Date("2000-01-01")
Locations = read.csv("../data/locations_area_20140807.csv")
Locations_original = Locations
House_area = mean(Locations_original$area[which(Locations_original$landuse=="HOUSE")])


#======================================================================================
# Parameters ----------
#======================================================================================
proto_pac = 5
# ratios of mosquito productivity per landuse relative to Houses (Morrisson 2006)
LandUseRatios = c(0.03, 0.01, 0.09, 1, 0.09, 0.09, 0.05, 0.003, 0.09,0.3) 

# Derived parameters
col.names = levels(Locations$landuse)
nlandusetypes = length(col.names)
colClasses = c(rep("numeric",nlandusetypes))

ProdR = read.table(text = "",
                    colClasses = colClasses,
                    col.names = col.names)  # table of landuseratios
ProdR[1,] = LandUseRatios  
Houses = which(Locations$landuse=="HOUSE")

# Including the Productivity ratio as an array in the Locations file 
ProdRatio = array(data = 0, dim = c(length(Locations[,1])))
for(i in 1:length(Locations[,1])){
  type = Locations$landuse[i]
  ProdRatio[i] = ProdR[,eval(type)] 
}
# Takes area into consideration by dividing the area by the avg house area
# Houses are assumed to have an average House area
# Other locations can't be more than 2 times a house

Locations$ProdRatioArea = ProdRatio * Locations$area /  House_area
Locations$ProdRatioArea[which(Locations$landuse == "HOUSE")] = 1
Locations$ProdRatioArea[which(Locations$ProdRatioArea > 2)] = 2

# Use any location to generate the trend:
SmallLocations = Locations[10000,]  


#======================================================================================
# Read temperature specific for Iquitos ------
#======================================================================================
Iquitos.climate = read.csv('../data/Iquitos_Climate_Bobby.csv')
startIndex = which(Iquitos.climate$date == as.character(tvec[1]))
# temperature = Iquitos.climate$temperature_mean[startIndex:(startIndex + length(dvec)-1) ]
# temperature.max = Iquitos.climate$temperature_max[startIndex:(startIndex + length(dvec)-1) ]
temperature = c()
temperature.max = c()
temperature.min = c()
for (i in 1:length(dvec)){
  temperature = c(temperature, rep(Iquitos.climate$temperature_mean[startIndex+i-1],length(tvec)/length(dvec)))
  temperature.max = c(temperature.max, rep(Iquitos.climate$temperature_max[startIndex+i-1],length(tvec)/length(dvec)))
  temperature.min = c(temperature.min, rep(Iquitos.climate$temperature_min[startIndex+i-1],length(tvec)/length(dvec)))
}
temp.mean = mean(temperature,na.rm = TRUE); temperature[is.na(temperature)] = temp.mean
temp.max.mean = mean(temperature.max,na.rm = TRUE); temperature.max[is.na(temperature.max)] = temp.max.mean
temp.min.mean = mean(temperature.min,na.rm = TRUE); temperature.min[is.na(temperature.min)] = temp.min.mean

sunexp = 0.1

water.temp.max = 15.03 + 0.27*temperature.min + 0.01*temperature.max^2 + 7.69*sunexp^2
water.temp.min = 5.02 - 1.36*sunexp + 0.81*temperature.min + 0.001*temperature.max^2
water.temp.mean = (water.temp.min + water.temp.max)/2
#======================================================================================
# obtain parameters governing death and development ------------
#======================================================================================
development_params = read.csv('../data/development_params.csv', sep = '\t', row.names=1)
mortality_thresholds = read.csv('../data/mortality_thresholds.csv', sep = '\t', row.names=1)

#Note that egg development is just used in the first time step, the full time series is estimated from the ODE system
egg_development_rate = calculate_development_rate_timeseries(development_params$eggs, water.temp.mean)*length(dvec)/length(tvec)
larval_development_rate = calculate_development_rate_timeseries(development_params$pupae, water.temp.mean)*length(dvec)/length(tvec)
pupal_development_rate = calculate_development_rate_timeseries(development_params$larvae, water.temp.mean)*length(dvec)/length(tvec)
gonotrophic_cycle_rate = calculate_development_rate_timeseries(development_params$adults, temperature)*length(dvec)/length(tvec)

egg_mortality = log(1/(1-calculate_mortality_rate_timeseries(mortality_thresholds$eggs, water.temp.max)*length(dvec)/length(tvec)))
#We need to somehow incorporate the density dependence of larval mortality, 
#but for now, use pupal rate for larvae (T dependent parts are identical) 
pupal_mortality = log(1/(1-calculate_mortality_rate_timeseries(mortality_thresholds$pupae, water.temp.max)*length(dvec)/length(tvec)))
larval_mortality = pupal_mortality
#DeathRate is adult mortality:
DeathRate = log(1/(1-calculate_adult_mortality_rate_timeseries(mortality_thresholds$adult, temperature.max)*length(dvec)/length(tvec)))
#old method:
#DeathRate = get_mortality_vector(tvec, temperature)

#======================================================================================
# Calculate Surface Space + Time for a small dataset --------
#======================================================================================
offsets = c("1.1","1.35","1.6","1.7","1.9", "2")
# small_abundance_surf = list()
o_f = offsets[2]
eval(parse(text = sprintf("load(\"EstimateSurface%s%s\")",o_f,".RData")))
XY = cbind(SmallLocations$xcor,SmallLocations$ycor)
daily_abundance_surf = Estimate.w.offset(XY,dvec,as.numeric(o_f))
normalized_daily_adults = exp(daily_abundance_surf[[1]])[1,]/exp(daily_abundance_surf[[1]])[1,1]
#Insert these into the midpoint of the day for the trend abundance surface
normalized_adults = rep(NA, length=length(tvec))
normalized_adults[seq(floor(length(tvec)/length(dvec)/2+1), length(tvec), length(tvec)/length(dvec))] <- normalized_daily_adults
adult_timeseries_trend = na.spline(normalized_adults)
# trend_abundance_surf = Estimate.w.offset.posix(XY,tvec,as.numeric(o_f)) # Log -- Est: [[1]] Std: [[2]]
# adult_timeseries_trend = exp(trend_abundance_surf[[1]])[1,]/exp(trend_abundance_surf[[1]])[1,1] #normalize by first timepoint to get non-location-specific trend
XYall = cbind(Locations$xcor, Locations$ycor)
log_initial_abundances = Estimate.w.offset(XYall,initial_date,as.numeric(o_f))
initial_abundances = proto_pac*exp(log_initial_abundances[[1]])[,1]
save(adult_timeseries_trend, initial_abundances, daily_abundance_surf, file="early_surfaces.RData")
#To get the trend for location i, do initial_abundances[i]*trend_abundance_surf

#==================================================
# Ways to deal with the nagative values
# 1. Make negative values 0 at end
# 2. Make negative values 0 throughout
# 3. Make negative values NA, fit a spline to the log, at end
# 4. Make negative values NA, fit a spline to the log, thoughout
# 5. Take a rolling average at end
# 6. Take a rolling average throughout
# 7. Do not make negative
#==================================================
#==============================
#The following tests different methods for calcuating the egg_dev elopment timeseries (no longer used):
#==============================
load("early_surfaces.RData")

larval_capacity = 0.2 #set to Inf for no density dependence
larval_power=1.
eggs_per_gon_cycle = 63

pupal_timeseries = get_pupal_timeseries(adult_timeseries_trend, DeathRate, pupal_development_rate, tvec)
save(pupal_timeseries, pupal_timeseries_nn2, pupal_timeseries_nn4, pupal_timeseries_nn6, file="pupal_timeseries.RData")


egg_timeseries = get_egg_timeseries_ode_2(adult_timeseries_trend,egg_mortality, egg_development_rate,gonotrophic_cycle_rate,eggs_per_gon_cycle, 15, tvec)
egg_timeseries_nn2 = pmax(egg_timeseries, rep(0, length(egg_timeseries)))
egg_timeseries_nn4 = exp(na.spline(log(egg_timeseries)))
egg_timeseries_nn6 = egg_timeseries
i=0
while (max(egg_timeseries_nn6<0)){
  i=i+1
  print(i)
  egg_timeseries_nn6 = rollapply(egg_timeseries, i, mean, partial = TRUE)
}
save(egg_timeseries, egg_timeseries_nn2, egg_timeseries_nn4, egg_timeseries_nn6, file="egg_timeseries.RData")


egg_timeseries = get_egg_timeseries_ode_2(adult_timeseries_trend,egg_mortality, egg_development_rate,gonotrophic_cycle_rate,eggs_per_gon_cycle, 15, tvec)
egg_timeseries_nn2 = pmax(egg_timeseries, rep(0, length(egg_timeseries)))
egg_timeseries_nn4 = exp(na.spline(log(egg_timeseries)))
egg_timeseries_nn6 = egg_timeseries
i=0
while (max(egg_timeseries_nn6<0)){
  i=i+1
  print(i)
  egg_timeseries_nn6 = rollapply(egg_timeseries, i, mean, partial = TRUE)
}
save(egg_timeseries, egg_timeseries_nn2, egg_timeseries_nn4, egg_timeseries_nn6, file="egg_timeseries.RData")

larval_timeseries = get_larval_timeseries_ode(egg_timeseries, pupal_timeseries, egg_development_rate, larval_development_rate, pupal_development_rate, larval_mortality, pupal_mortality, larval_capacity, larval_power, 1, tvec)
####UPTO HERE DONE OF NEW APPROACH#####


# Methods to compare========================================
# 1. Make negative values 0 at end
# 2. Make negative values 0 throughout
# 3. Make negative values NA, fit a spline to the log, at end
# 4. Make negative values NA, fit a spline to the log, thoughout
# 6. Take a rolling average throughout
# 7. Do not make negative
#===========================================================
population_timeseries = matrix(data=0, nrow=length(adult_timeseries_trend), ncol=7)
initial_conditions_1 = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
initial_conditions_2 = c(E =egg_timeseries_nn2[1], L = larval_timeseries_nn2[1], P = pupal_timeseries_nn2[1], N = adult_timeseries_trend[1])
initial_conditions_3 = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
initial_conditions_4 = c(E =egg_timeseries_nn4[1], L = larval_timeseries_nn4[1], P = pupal_timeseries_nn4[1], N = adult_timeseries_trend[1])
initial_conditions_5 = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
initial_conditions_6 = c(E =egg_timeseries_nn6[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])
initial_conditions_7 = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
params_1 = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn1, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_2 = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_3 = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn3, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_4 = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn4, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_5 = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn5, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_6 = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_7 = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
#population_timeseries[,1] = ode(y = initial_conditions_1, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_1)[,"N"]
#population_timeseries[,2] = ode(y = initial_conditions_2, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_2)[,"N"]
#population_timeseries[,3] = ode(y = initial_conditions_3, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_3)[,"N"]
#population_timeseries[,4] = ode(y = initial_conditions_4, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_4)[,"N"]
#population_timeseries[,5] = ode(y = initial_conditions_5, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_5)[,"N"]
population_timeseries[,6] = ode(y = initial_conditions_6, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_6)[,"N"]
#population_timeseries[,7] = ode(y = initial_conditions_7, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_7)[,"N"]
##6 is the best.

abs_errors = vector(length=7)
for (i in seq(length(abs_errors))){
  abs_errors[i] = sum(abs(adult_timeseries_trend-population_timeseries[,i])/adult_timeseries_trend)
}

#==================================================
# comparing with and without spraying -------------
#==================================================
spraying_thoroughness = 0.5
spray_today = vector(length = length(DeathRate))
start_day = 250; end_day= start_day+90
spray_today[start_day:end_day] = TRUE 

initial_conditions = c(E =egg_timeseries_nn6[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])
params_no_spraying = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_spraying = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate+spraying_thoroughness*spray_today, larval_capacity, larval_power)
population_timeseries_no_spraying = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_no_spraying)[,"N"]
population_timeseries_spraying = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_spraying)[,"N"]

plot(seq(1,length(tvec),1)[200:600], population_timeseries_spraying[200:600], type='l', ylim=c(0,1.5))
lines(seq(1,length(tvec),1)[200:600], population_timeseries_no_spraying[200:600], col='blue')
abline(v=start_day, col='green')
abline(v=end_day, col='red')

#==================================================
# from zero -------------------------
#==================================================
initial_conditions = c(E =0, L = 0, P = 0, N = 0.00000000001)
params = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params)[,"N"]

plot(seq(1,length(tvec),1)[1:365], population_timeseries[1:365], type='l')
lines(adult_timeseries_trend, col='blue')

#==================================================
# check the un-normalized equations
#==================================================
locations_to_test = floor(runif(3,0,length(initial_abundances)))+1
initial_conditions = list()
params = initial_conditions
population_timeseries = matrix(data=0, nrow=length(adult_timeseries_trend), ncol=length(locations_to_test))
for (i in seq(length(locations_to_test))){
  initial_conditions[[i]] = c(E =egg_timeseries_nn6[1]*initial_abundances[locations_to_test[i]], L = larval_timeseries_nn6[1]*initial_abundances[locations_to_test[i]], P = pupal_timeseries_nn6[1]*initial_abundances[locations_to_test[i]], N = adult_timeseries_trend[1]*initial_abundances[locations_to_test[i]])
  params[[i]] = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity*initial_abundances[locations_to_test[i]]^larval_power, larval_power)
  population_timeseries[,i] = ode(y = initial_conditions[[i]], times = seq(1,length(tvec),1), func = full_model_ode, parms = params[[i]])[,"N"]
}

plot(population_timeseries[,1]/population_timeseries[,2], type='l')

spraying_thoroughness = 0.5
spray_today = vector(length = length(DeathRate))
start_day = 250; end_day= start_day+90
spray_today[start_day:end_day] = TRUE 
initial_conditions = list()
params_spraying = initial_conditions
population_timeseries_spraying = matrix(data=0, nrow=length(adult_timeseries_trend), ncol=length(locations_to_test))
params_no_spraying = initial_conditions
population_timeseries_no_spraying = matrix(data=0, nrow=length(adult_timeseries_trend), ncol=length(locations_to_test))
for (i in seq(length(locations_to_test))){
  initial_conditions[[i]] = c(E =egg_timeseries_nn6[1]*initial_abundances[locations_to_test[i]], L = larval_timeseries_nn6[1]*initial_abundances[locations_to_test[i]], P = pupal_timeseries_nn6[1]*initial_abundances[locations_to_test[i]], N = adult_timeseries_trend[1]*initial_abundances[locations_to_test[i]])
  params_no_spraying[[i]] = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity*initial_abundances[locations_to_test[i]]^larval_power, larval_power)
  params_spraying[[i]] = list(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate+spraying_thoroughness*spray_today, larval_capacity*initial_abundances[locations_to_test[i]]^larval_power, larval_power)
  population_timeseries_no_spraying[,i] = ode(y = initial_conditions[[i]], times = seq(1,length(tvec),1), func = full_model_ode, parms = params_no_spraying[[i]])[,"N"]
  population_timeseries_spraying[,i] = ode(y = initial_conditions[[i]], times = seq(1,length(tvec),1), func = full_model_ode, parms = params_spraying[[i]])[,"N"]
}
png("/home/sean/Desktop/perturbed_ratio.png")
plot(population_timeseries_spraying[,1]/population_timeseries_spraying[,2], type='l')
dev.off()
plot(population_timeseries_spraying[,1], typ='l')
lines(population_timeseries_spraying[,2])
abline(v=start_day, col='green')
abline(v=end_day, col='red')



#======================================================================================
# Write to output---------
#======================================================================================
# First, the location varying parameters:
Locations$initial_adults = initial_abundances
save_columns = c("xcor","ycor","landuse","block","zone","neighborhood","code","initial_adults")
LocationsMosquitoes = Locations[save_columns]
write.csv(LocationsMosquitoes,"../output/Locations20160720.csv",row.names = F, quote = FALSE)

# save(Emergence_TS,file = sprintf("../data/Emergence%s.RData",o_f))

# Now the fixed parameters:
fileName <- file("../output/fixed_parameters.txt")
open(fileName, "wt")
writeLines(paste("initial_pupae_normalized = ", pupal_timeseries_nn6[1]), fileName)
writeLines(paste("initial_larvae_normalized = ", larval_timeseries_nn6[1]), fileName)
writeLines(paste("initial_eggs_normalized = ", egg_timeseries_nn6[1]), fileName)
writeLines(paste("larval_capacity_normalized = ", larval_capacity), fileName)
writeLines(paste("larval_power = ", larval_power), fileName)
writeLines(paste("total_initial_population = ", sum(initial_abundances)), fileName)
close(fileName)

# The time-varying parameters are done in generate_parameters_climate. Save them here:
save(egg_development_rate, larval_development_rate, pupal_development_rate, female_egg_laying_timeseries_nn6, egg_mortality, larval_mortality, pupal_mortality, DeathRate, file="time_dependent_series.RData")

