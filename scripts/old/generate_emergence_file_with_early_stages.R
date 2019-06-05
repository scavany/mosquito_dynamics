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

pupal_mortality = log(1/(1-calculate_mortality_rate_timeseries(mortality_thresholds$pupae, water.temp.max)*length(dvec)/length(tvec)))
larval_mortality = pupal_mortality

larval_capacity = 0.2 #set to Inf for no density dependence
larval_power=1.
eggs_per_gon_cycle = 63
correction_factor = 0.985

#additional_random_death_rate = 0.0
#pupal_mortality = log(1/(1-calculate_mortality_rate_timeseries(mortality_thresholds$pupae, water.temp.max)*length(dvec)/length(tvec)))
#pupal_mortality = pupal_mortality + additional_random_death_rate
#larval_mortality = pupal_mortality


pupal_timeseries = get_pupal_timeseries(adult_timeseries_trend, DeathRate, pupal_development_rate, tvec)
pupal_timeseries_nn2 = pmax(pupal_timeseries, rep(0, length(pupal_timeseries)))
pupal_timeseries_nn4 = exp(na.spline(log(pupal_timeseries)))
pupal_timeseries_nn6 = pupal_timeseries
i=0
while (max(pupal_timeseries_nn6<0)){
  i=i+1
  print(i)
  pupal_timeseries_nn6 = rollapply(pupal_timeseries, i, mean, partial = TRUE)
}
pupal_test_pop = ode(y = c(N=1), times = seq(1,length(tvec),1), func = pupal_test_ode, parms = list(pupal_development_rate, pupal_timeseries, DeathRate), method='rk')[,"N"]
save(pupal_timeseries, pupal_timeseries_nn2, pupal_timeseries_nn4, pupal_timeseries_nn6, file="pupal_timeseries.RData")



larval_timeseries = get_larval_timeseries(pupal_timeseries, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_timeseries_nn2 = get_larval_timeseries(pupal_timeseries_nn2, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_timeseries_nn4 = get_larval_timeseries(pupal_timeseries_nn4, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_timeseries_nn6_initial = get_larval_timeseries(pupal_timeseries_nn6, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_timeseries_nn2 = pmax(larval_timeseries_nn2, rep(0, length(larval_timeseries_nn2)))
larval_timeseries_nn4 = exp(na.spline(log(larval_timeseries_nn4)))
larval_timeseries_nn6 = larval_timeseries_nn6_initial
i=0
while (max(larval_timeseries_nn6<0)){
  i=i+1
  print(i)
  larval_timeseries_nn6 = rollapply(larval_timeseries_nn6_initial, i, mean, partial = TRUE)
}
larval_test_pop = ode(y = c(P=pupal_timeseries[1], N=adult_timeseries_trend[1]), times = seq(1,length(tvec),1), func = larval_test_ode, parms = list(larval_development_rate, pupal_development_rate, larval_timeseries, pupal_mortality, DeathRate))[,"N"]
save(larval_timeseries, larval_timeseries_nn2, larval_timeseries_nn4, larval_timeseries_nn6, file="larval_timeseries.RData")

#The next step needs an improved larval mortality (incorporating density dependence)
larval_inflow_timeseries = get_larval_inflow_timeseries(larval_timeseries, larval_mortality, larval_development_rate, larval_capacity, larval_power, tvec)
larval_inflow_timeseries_nn2 = get_larval_inflow_timeseries(larval_timeseries_nn2, larval_mortality, larval_development_rate, larval_capacity, larval_power, tvec)
larval_inflow_timeseries_nn4 = get_larval_inflow_timeseries(larval_timeseries_nn4, larval_mortality, larval_development_rate, larval_capacity, larval_power, tvec)
larval_inflow_timeseries_nn6_initial = get_larval_inflow_timeseries(larval_timeseries_nn6, larval_mortality, larval_development_rate, larval_capacity, larval_power, tvec)
larval_inflow_timeseries_nn2 = pmax(larval_inflow_timeseries_nn2, rep(0, length(larval_inflow_timeseries_nn2)))
larval_inflow_timeseries_nn4 = exp(na.spline(log(larval_inflow_timeseries_nn4)))
larval_inflow_timeseries_nn6 = larval_inflow_timeseries_nn6_initial
i=0
while (max(larval_inflow_timeseries_nn6<0)){
  i=i+1
  print(i)
  larval_inflow_timeseries_nn6 = rollapply(larval_inflow_timeseries_nn6_initial, i, mean, partial = TRUE)
}
larval_inflow_test_pop = ode(y = c(L= larval_timeseries[1], P=pupal_timeseries[1], N=adult_timeseries_trend[1]), times = seq(1,length(tvec),1), func = larval_inflow_test_ode, parms = list(larval_development_rate, pupal_development_rate, larval_inflow_timeseries, larval_mortality,  pupal_mortality, DeathRate, larval_capacity, larval_power))[,"N"]
save(larval_inflow_timeseries, larval_inflow_timeseries_nn2, larval_inflow_timeseries_nn4, larval_inflow_timeseries_nn6, file="larval_inflow_timeseries_Inf_1.RData")

##NOTE: All forecasts seem to be good using timeseries up until this point.

# #Get egg timeseries
# egg_timeseries = get_egg_timeseries(larval_timeseries, larval_mortality, larval_development_rate, egg_development_rate, larval_capacity, larval_power, tvec)
# egg_timeseries_nn2 = get_egg_timeseries(larval_timeseries_nn2, larval_mortality, larval_development_rate, egg_development_rate, larval_capacity, larval_power, tvec)
# egg_timeseries_nn4 = get_egg_timeseries(larval_timeseries_nn4, larval_mortality, larval_development_rate, egg_development_rate, larval_capacity, larval_power, tvec)
# egg_timeseries_nn6_initial = get_egg_timeseries(larval_timeseries_nn6, larval_mortality, larval_development_rate, egg_development_rate, larval_capacity, larval_power, tvec)
# egg_timeseries_nn2 = pmax(egg_timeseries_nn2, rep(0, length(egg_timeseries_nn2)))
# egg_timeseries_nn4 = exp(na.spline(log(egg_timeseries_nn4)))
# egg_timeseries_nn6 = egg_timeseries_nn6_initial
# i=0
# while (max(egg_timeseries_nn6<0)){
#   i=i+1
#   print(i)
#   egg_timeseries_nn6 = rollapply(egg_timeseries_nn6_initial, i*24, mean, partial = TRUE)
# }
# save(egg_timeseries, egg_timeseries_nn2, egg_timeseries_nn4, egg_timeseries_nn6, file="egg_timeseries_Inf_1.RData")


# #Get gonotrophcycle rate
# #This is the timeseries we will use as input to the ABM:
# female_egg_laying_timeseries = get_female_egg_laying_rate(adult_surf.in = adult_timeseries_trend,egg_surf.in = egg_timeseries,egg_mortality.in =egg_mortality ,egg_dev.in = egg_development_rate,times.in = tvec)
# female_egg_laying_timeseries_nn1 = pmax(female_egg_laying_timeseries, rep(0, length(female_egg_laying_timeseries)))
# female_egg_laying_timeseries_nn2 = get_female_egg_laying_rate(adult_surf.in = adult_timeseries_trend,egg_surf.in = egg_timeseries_nn2,egg_mortality.in =egg_mortality ,egg_dev.in = egg_development_rate,times.in = tvec)
# female_egg_laying_timeseries_nn2 = pmax(female_egg_laying_timeseries_nn2, rep(0, length(female_egg_laying_timeseries_nn2)))
# female_egg_laying_timeseries_nn3 = exp(na.spline(log(female_egg_laying_timeseries)))
# female_egg_laying_timeseries_nn4 = get_female_egg_laying_rate(adult_surf.in = adult_timeseries_trend,egg_surf.in = egg_timeseries_nn4,egg_mortality.in =egg_mortality ,egg_dev.in = egg_development_rate,times.in = tvec)
# female_egg_laying_timeseries_nn4 = exp(na.spline(log(female_egg_laying_timeseries_nn4)))
# female_egg_laying_timeseries_nn5 = female_egg_laying_timeseries
# i=0
# while (max(female_egg_laying_timeseries_nn5<0)){
#   i=i+1
#   print(i)
#   female_egg_laying_timeseries_nn5 = rollapply(female_egg_laying_timeseries, i, mean, partial = TRUE)
# }
# female_egg_laying_timeseries_nn6_initial = get_female_egg_laying_rate(adult_surf.in = adult_timeseries_trend,egg_surf.in = egg_timeseries_nn6,egg_mortality.in =egg_mortality ,egg_dev.in = egg_development_rate,times.in = tvec)
# female_egg_laying_timeseries_nn6 = female_egg_laying_timeseries_nn6_initial
# i=0
# while (max(female_egg_laying_timeseries_nn6<0)){
#   i=i+1
#   print(i)
#   female_egg_laying_timeseries_nn6 = rollapply(female_egg_laying_timeseries_nn6_initial, i, mean, partial = TRUE)
# }
# save(female_egg_laying_timeseries, female_egg_laying_timeseries_nn1, female_egg_laying_timeseries_nn2, female_egg_laying_timeseries_nn3, female_egg_laying_timeseries_nn4, female_egg_laying_timeseries_nn5, female_egg_laying_timeseries_nn6, file="female_egg_laying_timeseries_Inf_1.RData")

eggs0 = larval_inflow_timeseries[1]/egg_development_rate[1]
eggs0_nn2 = max(0,larval_inflow_timeseries_nn2[1]/egg_development_rate[1])
eggs0_nn4 = max(0,larval_inflow_timeseries_nn4[1]/egg_development_rate[1])
eggs0_nn6 = max(0,larval_inflow_timeseries_nn6[1]/egg_development_rate[1])

egg_timeseries = get_egg_timeseries_ode(larval_inflow_timeseries, adult_timeseries_trend, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0, tvec)
egg_timeseries_nn2 = get_egg_timeseries_ode(larval_inflow_timeseries_nn2, adult_timeseries_trend, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0_nn2, tvec)
egg_timeseries_nn4 = get_egg_timeseries_ode(larval_inflow_timeseries_nn4, adult_timeseries_trend, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0_nn4, tvec)
egg_timeseries_nn6_initial = get_egg_timeseries_ode(larval_inflow_timeseries_nn6, adult_timeseries_trend, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0_nn6, tvec)
egg_timeseries_nn2 = pmax(egg_timeseries_nn2, rep(0, length(egg_timeseries_nn2)))
egg_timeseries_nn4 = exp(na.spline(log(egg_timeseries_nn4)))
egg_timeseries_nn6 = egg_timeseries_nn6_initial
i=0
while (max(egg_timeseries_nn6<0)){
  i=i+1
  print(i)
  egg_timeseries_nn6 = rollapply(egg_timeseries_nn6_initial, i, mean, partial = TRUE)
}


egg_development_timeseries = larval_inflow_timeseries/egg_timeseries
egg_development_timeseries_nn1 = pmax(egg_development_timeseries, rep(0, length(egg_development_timeseries)))
egg_development_timeseries_nn2 = larval_inflow_timeseries_nn2/egg_timeseries_nn2
egg_development_timeseries_nn3 = exp(na.spline(log(egg_development_timeseries)))
egg_development_timeseries_nn4 = larval_inflow_timeseries_nn4/egg_timeseries_nn4
#egg_development_timeseries_nn5 = egg_development_timeseries
#i=0
#while (max(egg_development_timeseries_nn5<0)){
#  i=i+1
#  print(i)
#  egg_development_timeseries_nn5 = rollapply(egg_development_timeseries, i, mean, partial = TRUE)
#}
egg_development_timeseries_nn6 = larval_inflow_timeseries_nn6/egg_timeseries_nn6                                                                                                                                                                                                         



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
#initial_conditions_5 = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
initial_conditions_6 = c(E =egg_timeseries_nn6[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])
initial_conditions_7 = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
params_1 = list(egg_development_timeseries_nn1, larval_development_rate, pupal_development_rate, correction_factor*gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_2 = list(egg_development_timeseries_nn2, larval_development_rate, pupal_development_rate, correction_factor*gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_3 = list(egg_development_timeseries_nn3, larval_development_rate, pupal_development_rate, correction_factor*gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_4 = list(egg_development_timeseries_nn4, larval_development_rate, pupal_development_rate, correction_factor*gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
#params_5 = list(egg_development_timeseries_nn5, larval_development_rate, pupal_development_rate, correction_factor*gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_6 = list(egg_development_timeseries_nn6, larval_development_rate, pupal_development_rate, correction_factor*gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
params_7 = list(egg_development_timeseries, larval_development_rate, pupal_development_rate, correction_factor*gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
population_timeseries[,1] = ode(y = initial_conditions_1, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_1, method='rk')[,"N"]
population_timeseries[,2] = ode(y = initial_conditions_2, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_2, method='rk')[,"N"]
population_timeseries[,3] = ode(y = initial_conditions_3, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_3, method='rk')[,"N"]
population_timeseries[,4] = ode(y = initial_conditions_4, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_4, method='rk')[,"N"]
#population_timeseries[,5] = ode(y = initial_conditions_5, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_5, method='rk')[,"N"]
population_timeseries[,6] = ode(y = initial_conditions_6, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_6, method='rk')[,"N"]
population_timeseries[,7] = ode(y = initial_conditions_7, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_7, method='rk')[,"N"]
##6 is the best.

abs_errors = vector(length=7)
for (i in seq(length(abs_errors))){
  abs_errors[i] = sum(abs(adult_timeseries_trend-population_timeseries[,i])/adult_timeseries_trend)
}

save(population_timeseries_1, population_timeseries_2, population_timeseries_3, population_timeseries_4, population_timeseries_5, population_timeseries_6, population_timeseries_7, file='egg_dev_comparison_Inf_1.RData')


#==================================================================
# Now, using the seventh method (ignoring negative values), which seems most accurate, see which larval capacity seems most accurate
# Also, are doing this with the normalized equation set - so location is not important
#==================================================================
# load("surfaces.RData")
# egg_timeseries_list=list()
# mean_errors = c()
# larval_cap = c()
# population_timeseries_list = list()
# 
# pupal_timeseries = get_pupal_timeseries(adult_timeseries_trend, DeathRate, pupal_development_rate, tvec)
# larval_timeseries = get_larval_timeseries(pupal_timeseries, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
# eggs_per_gon_cycle = 63 #Otero, 2006
# larval_power=1.
# # In the below loop, values of i below -3 fail
# j=1
# for (larval_capacity in c(10^(-5:5), Inf)){  
#   print(j)
#   #The next step needs an improved larval mortality (incorporating density dependence)
#   larval_inflow_timeseries = get_larval_inflow_timeseries(larval_timeseries, larval_mortality, larval_development_rate, larval_capacity, larval_power, tvec)
#   
#   #estimate the initial eggs, and then the egg timeseries:
# 
#   eggs0 = larval_inflow_timeseries[1]/egg_development_rate[1]
#   egg_timeseries = get_egg_timeseries_ode(larval_inflow_timeseries, adult_timeseries_trend, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0, tvec)
#   # 
#   # #This is the timeseries we will use as input to the ABM:
#   egg_development_timeseries = larval_inflow_timeseries/egg_timeseries
#   
#   initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
#   params = list(egg_development_timeseries, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power) 
#   population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params)
#   
#   egg_timeseries_list[[j]] = egg_development_timeseries
#   population_timeseries_list[[j]] = population_timeseries
#   mean_errors = append(mean_errors, mean(abs((adult_timeseries_trend-population_timeseries[,"N"])/adult_timeseries_trend)))
#   larval_cap = append(larval_cap, larval_capacity)
#   mean(abs((adult_timeseries_trend-population_timeseries[,"N"])/adult_timeseries_trend))
#   # plot(population_timeseries)
#   # dev.copy(png, paste0("../figures/data_timeseries_",larval_capacity, ".png"))
#   # dev.off()
#   
#   j=j+1
# }
# # Best order of magnitude is around 1, with mean errors of 0.3641537
# # Best is with larval_capacity of around 0.7, with mean errors 0.1788495
# 
# save(population_timeseries_list, mean_errors, egg_timeseries_list, larval_cap, file="larval_capacity_comparison_order.Rdata")

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
writeLines(paste("initial_pupae_normalized = ", pupal_timeseries[1]), fileName)
writeLines(paste("initial_larvae_normalized = ", larval_timeseries[1]), fileName)
writeLines(paste("initial_eggs_normalized = ", egg_timeseries[1]), fileName)
writeLines(paste("larval_capacity_normalized = ", larval_capacity), fileName)
writeLines(paste("larval_power = ", larval_power), fileName)
writeLines(paste("eggs_per_gon_cycle = ", eggs_per_gon_cycle*correction_factor), fileName)
close(fileName)


fixed_parameters = data.frame("initial_pupae_normalized"=c(pupal_timeseries[1])) 
fixed_parameters$initial_larvae_normalized = c(larval_timeseries[1])
fixed_parameters$initial_eggs_normalized = c(egg_timeseries[1])
fixed_parameters$larval_capacity_normalized = c(larval_capacity)
fixed_parameters$larval_power = c(larval_power)
fixed_parameters$eggs_per_gon_cycle = c(eggs_per_gon_cycle*correction_factor)
write.csv(fixed_parameters, "../output/fixed_parameters.csv", row.names=F, quote = FALSE)

# The time-varying parameters are done in generate_parameters_climate. Save them here:
save(egg_development_timeseries_nn6, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, DeathRate, file="time_dependent_series.RData")

