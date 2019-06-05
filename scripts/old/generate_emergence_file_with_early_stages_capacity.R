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

source('./generate_emergence_functions_with_early_stages_capacity.R')
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
save(adult_timeseries_trend, initial_abundances, daily_abundance_surf, file="surfaces_lc.RData")
#To get the trend for location i, do initial_abundances[i]*trend_abundance_surf

load("surfaces_lc_extended.RData")

#======================================================================================
# GEt timeseries
#======================================================================================
larval_power=1.0
eggs_per_gon_cycle = 63
adult_development_rate = eggs_per_gon_cycle*gonotrophic_cycle_rate/2
correction_factor = 0.242
eggs0 = 18.4#this is approximately right, by inspection

pupal_timeseries = get_pupal_timeseries(adult_timeseries_trend, DeathRate, pupal_development_rate, tvec)
adult_ts_check_pupae = ode(y=adult_timeseries_trend[1], times=seq(1,length(tvec),1), func = pupal_test_ode, parms = list(pupal_timeseries, DeathRate, pupal_development_rate), method='rk')

larval_timeseries = get_larval_timeseries(pupal_timeseries, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
pupal_ts_check_larvae = ode(y=pupal_timeseries[1], times=seq(1,length(tvec),1), func = larval_test_ode, parms = list(larval_timeseries, pupal_mortality, pupal_development_rate, larval_development_rate), method='rk')
adult_ts_check_larvae = ode(y=adult_timeseries_trend[1], times=seq(1,length(tvec),1), func = pupal_test_ode, parms = list(pupal_ts_check_larvae[,2], DeathRate, pupal_development_rate), method='rk')

egg_timeseries = get_egg_timeseries(adult_timeseries_trend, egg_mortality, adult_development_rate, egg_development_rate, eggs0, tvec)
adult_ts_check_egg = (differentiate_timeseries(egg_timeseries, tvec) + egg_mortality*egg_timeseries + egg_development_rate*egg_timeseries)/adult_development_rate

larval_capacity_timeseries = get_larval_capacity_timeseries(larval_timeseries, egg_timeseries, larval_mortality, larval_development_rate, egg_development_rate,larval_power, tvec)

initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
params = list(DeathRate, pupal_mortality, larval_mortality, egg_mortality, correction_factor*adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity_timeseries, larval_power)
population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params)
#I tried different larval_powers - it does not seem to make much difference. for simiplicity, stick with larval_power =1

minimize_function = function(correction_factor, DeathRate, pupal_mortality, larval_mortality, egg_mortality, adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity_timeseries, larval_power){
  params = list(DeathRate, pupal_mortality, larval_mortality, egg_mortality, correction_factor*adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity_timeseries, larval_power)
  population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params)
  return(abs(mean(population_timeseries[,"N"]) - mean(adult_timeseries_trend)))
}

correction_factor = optimize(minimize_function, DeathRate, pupal_mortality, larval_mortality, egg_mortality, adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity_timeseries, larval_power, interval=c(0,1))
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
fileName <- file("../output/fixed_parameters_capacity.txt")
open(fileName, "wt")
writeLines(paste("initial_pupae_normalized = ", pupal_timeseries[1]), fileName)
writeLines(paste("initial_larvae_normalized = ", larval_timeseries[1]), fileName)
writeLines(paste("initial_eggs_normalized = ", egg_timeseries[1]), fileName)
writeLines(paste("larval_power = ", larval_power), fileName)
writeLines(paste("eggs_per_gon_cycle = ", eggs_per_gon_cycle*correction_factor$minimum), fileName)
close(fileName)


# The time-varying parameters are done in generate_parameters_climate. Save them here:
save(larval_capacity_timeseries, egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, DeathRate, file="time_dependent_series_capacity.RData")

