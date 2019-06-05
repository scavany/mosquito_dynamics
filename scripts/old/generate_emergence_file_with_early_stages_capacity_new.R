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
o_f = offsets[6]
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
additional_random_death_rate = log(1/(1-0.11/3)) #1/3 of containers randomly filled, 0.11 probability per day of being emptied
# additional_random_death_rate = 0.352 # This is chosen so that the gonotrophic cycle rate matches that given by skeeterbuster and assuming 63 eggs per gc
pupal_mortality = log(1/(1-calculate_mortality_rate_timeseries(mortality_thresholds$pupae, water.temp.max)*length(dvec)/length(tvec)))
pupal_mortality = pupal_mortality + additional_random_death_rate
larval_mortality = pupal_mortality

larval_power = 1
eggs_per_gon_cycle = 63

pupal_timeseries = (differentiate_timeseries(adult_timeseries_trend, tvec) + DeathRate*adult_timeseries_trend)/pupal_development_rate

egg_timeseries = get_egg_timeseries_ode_2(adult_timeseries_trend, egg_mortality, egg_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, 15, tvec, method='euler')

larval_timeseries = (differentiate_timeseries(pupal_timeseries, tvec) + pupal_timeseries*(pupal_mortality + pupal_development_rate))/larval_development_rate

larval_capacity_timeseries = larval_timeseries^larval_power*(differentiate_timeseries(larval_timeseries, tvec) - egg_development_rate*egg_timeseries + larval_timeseries*(larval_mortality + larval_development_rate))^-1


extra_mortality_nn6 = (egg_development_rate*egg_timeseries - differentiate_timeseries(larval_timeseries_nn6, tvec))/larval_timeseries_nn6 - larval_mortality - larval_development_rate - larval_timeseries_nn6^larval_power/larval_capacity
extra_mortality_nn6_initial = extra_mortality_nn6 
i=0
while (max(extra_mortality_nn6 <0)){
  i=i+1
  print(i)
  extra_mortality_nn6 = rollapply(extra_mortality_nn6_initial, i, mean, partial = TRUE)
}


initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries_trend[1])
params = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate*eggs_per_gon_cycle/2, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity, larval_power)
population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params, method='euler')[,"N"]


###correction factor is not necessaryily needed - but see what fits best anyways
initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])

minimize_function = function(correction_factor, DeathRate, pupal_mortality, larval_mortality, egg_mortality, extra_mortality_nn6, adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity, larval_power, eggs_per_gon_cycle){
  params = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, correction_factor*eggs_per_gon_cycle)
  population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params)
  return(sum(abs(adult_timeseries_trend-population_timeseries[,"N"])/adult_timeseries_trend))
}

correction_factor_optimized = optimize(minimize_function, DeathRate, pupal_mortality, larval_mortality, egg_mortality, extra_mortality_nn6, adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity, larval_power, eggs_per_gon_cycle, interval=c(0.8,1.25))

correction_factor = correction_factor_optimized$minimum
params = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, correction_factor*eggs_per_gon_cycle)
population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params)
plot(population_timeseries[,"N"])
lines(adult_timeseries_trend, col='red')
#correction factor is 1.05 - probably can do without it?

#==================================================
# comparing with and without spraying -------------
#==================================================
spraying_thoroughness = 0.5
spray_today = vector(length = length(DeathRate))
start_day = 250; end_day= start_day+90
spray_today[start_day:end_day] = TRUE 

initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])
params_no_spraying = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, correction_factor_optimized$minimum*eggs_per_gon_cycle)
params_spraying = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate+spraying_thoroughness*spray_today, larval_capacity, larval_power, correction_factor_optimized$minimum*eggs_per_gon_cycle)
population_timeseries_no_spraying = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_no_spraying)[,"N"]
population_timeseries_spraying = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_spraying)[,"N"]

plot(seq(1,length(tvec),1)[200:600], population_timeseries_spraying[200:600], type='l', ylim=c(0,1.5))
lines(seq(1,length(tvec),1)[200:600], population_timeseries_no_spraying[200:600], col='blue')
abline(v=start_day, col='green')
abline(v=end_day, col='red')

#==================================================
# trying a prediction with log equations
#==================================================
egg_ts = vector(length=length(adult_timeseries_trend))
larval_ts = vector(length=length(adult_timeseries_trend))
pupal_ts = vector(length=length(adult_timeseries_trend))
adult_ts = vector(length=length(adult_timeseries_trend))
egg_ts[1] = egg_timeseries[1]
larval_ts[1] = larval_timeseries_nn6[1]
pupal_ts[1] = pupal_timeseries_nn6[1]
adult_ts[1] = adult_timeseries_trend[1]

for (i in seq(2, length(adult_timeseries_trend))){
  print(i)
  egg_ts[i] = exp(log(egg_ts[i-1]) + adult_ts[i-1]*gonotrophic_cycle_rate[i-1]*eggs_per_gon_cycle/egg_ts[i-1]/2 - egg_mortality[i-1] - egg_development_rate[i-1])
  larval_ts[i] = exp(log(larval_ts[i-1]) + egg_ts[i-1]*egg_development_rate[i-1]/larval_ts[i-1] - larval_mortality[i-1] - larval_development_rate[i-1] - extra_mortality_nn6[i-1] - larval_ts[i-1]^larval_power/larval_capacity)
  pupal_ts[i] = exp(log(pupal_ts[i-1]) + larval_development_rate[i-1]*larval_ts[i-1] - pupal_mortality[i-1] - pupal_development_rate[i-1] - extra_mortality_nn6[i-1])
  adult_ts[i] = exp(log(adult_ts[i-1]) + pupal_development_rate[i-1]*pupal_ts[i-1] - DeathRate[i-1])
}


#==================================================
# from zero -------------------------
#==================================================
initial_conditions = c(E =1e-10, L = 0, P = 0, N = 0)
params = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, correction_factor*eggs_per_gon_cycle)
population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params)[,"N"]

plot(seq(1,length(tvec),1)[1:365], population_timeseries[1:365], type='l')
lines(adult_timeseries_trend, col='blue')

#==================================================
# does the density dependence maintain spatial heterogeneity? 
#=================================================
#Look at the extremes, start with 0 adults
maxN0 = max(initial_abundances)
minN0 = min(initial_abundances)

initial_conditions = c(E =1e-30, L = 0, P = 0, N = 0)
params_max = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity*maxN0^larval_power, larval_power, correction_factor*eggs_per_gon_cycle)
params_min = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity*minN0^larval_power, larval_power, correction_factor*eggs_per_gon_cycle)
population_timeseries_max = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_max)[,"N"]
population_timeseries_min = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_min)[,"N"]

png("../figures/larval_capacity_equilibrium.png")
plot(tvec[1:1000], population_timeseries_max[1:1000], type='l', ylim = c(0,8))
lines(tvec[1:1000], adult_timeseries_trend[1:1000]*maxN0, col='black', lty='dashed')
lines(tvec[1:1000], adult_timeseries_trend[1:1000]*minN0, col='red', lty='dashed')
lines(tvec[1:1000], population_timeseries_min[1:1000], col='red')
legend(x=as.Date("2000-01-01"), y=8, legend=c("largest location, prediction","largest location, data","smallest location, prediction","smallest location, data"), col=c("black", "black", "red", "red"), lty = c("solid", "dashed", "solid", "dashed"))
dev.off()


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
writeLines(paste("initial_eggs_normalized = ", egg_timeseries[1]), fileName)
writeLines(paste("larval_capacity_normalized = ", larval_capacity), fileName)
writeLines(paste("larval_power = ", larval_power), fileName)
#writeLines(paste("total_initial_population = ", sum(initial_abundances)), fileName)
close(fileName)

# The time-varying parameters are done in generate_parameters_climate. Save them here:
save(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, DeathRate, extra_mortality_nn6, correction_factor, eggs_per_gon_cycle, file="time_dependent_series.RData")

