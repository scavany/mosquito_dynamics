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
require(sdprisk)
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
#gonotrophic_cycle_rate = calculate_development_rate_timeseries(development_params$adults, temperature)*length(dvec)/length(tvec)

#egg_mortality = calculate_mortality_rate_timeseries(mortality_thresholds$eggs, water.temp.max)*length(dvec)/length(tvec)
#We need to somehow incorporate the density dependence of larval mortality, 
#but for now, use pupal rate for larvae (T dependent parts are identical) 
#pupal_mortality = calculate_mortality_rate_timeseries(mortality_thresholds$pupae, water.temp.max)*length(dvec)/length(tvec)
#larval_mortality = pupal_mortality
#DeathRate is adult mortality:
DeathRate = calculate_adult_mortality_rate_timeseries(mortality_thresholds$adult, temperature.max)*length(dvec)/length(tvec)
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
save(adult_timeseries_trend, initial_abundances, daily_abundance_surf, file="surfaces.RData")
#To get the trend for location i, do initial_abundances[i]*trend_abundance_surf

#======================================================================================
# Get pdfs etc.
#======================================================================================
load("surfaces.RData")
times=seq(length(tvec))
times0=seq(length(tvec)-1)

Lambda_E = Lambda(times, times0, egg_development_rate)
Lambda_L = Lambda(times, times0, larval_development_rate)
Lambda_P = Lambda(times, times0, pupal_development_rate)

f_E = dnonhomexp(times, times0, egg_development_rate, Lambda_E)   
f_L = dnonhomexp(times, times0, larval_development_rate, Lambda_L)   
f_P = dnonhomexp(times, times0, pupal_development_rate, Lambda_P)   

f_M = convolution(f_E, f_L)
f_N = convolution(f_M, f_P)

save(f_E, f_L, f_P, f_M, f_N, file = "non_homogeneous_exponential_densities.RData")
#======================================================================================
# Get timeseries
#======================================================================================
lambda_e = mean(egg_development_rate)
lambda_l = mean(larval_development_rate)
lambda_p = mean(pupal_development_rate)
lambdas = c(lambda_e, lambda_l, lambda_p)

dNdt = differentiate_timeseries(adult_timeseries_trend, tvec)
#almost all the eggs will have matured by:
maxn = floor(qhypoexp(1-1e-8, lambdas))
times = seq(length(tvec))

#do the integration, using trapezium rule
number_of_strips = maxn*5
density_integral_mean = integrate_trapezoid(integrand=dhypoexp, lower=0, upper=maxn, nstrips=number_of_strips, rate=lambdas)
print(density_integral_mean) #Should be 1 - using strips 5x maxn gives a density of 1
density_integral_varying = integrate_trapezoid(integrand=density_integrand_varying, lower=0, upper=maxn, nstrips=number_of_strips, times, egg_development_rate, larval_development_rate, pupal_development_rate)

emergence_integral_mean = integrate_trapezoid(integrand=emergence_integrand_mean, lower=0, upper=maxn, nstrips = number_of_strips, times, adult_timeseries_trend, lambdas )
emergence_integral_varying = integrate_trapezoid(integrand=emergence_integrand_varying, lower=0, upper=maxn, nstrips = number_of_strips, times, adult_timeseries_trend, egg_development_rate, larval_development_rate, pupal_development_rate)

emergence_mean = (dNdt + DeathRate*adult_timeseries_trend)*density_integral_mean/emergence_integral_mean
emergence_varying = (dNdt + DeathRate*adult_timeseries_trend)*density_integral_varying/emergence_integral_varying

save(density_integral_varying, emergence_integral_mean, emergence_integral_varying, emergence_mean, emergence_varying, file="delay_integrals.RData")

#======================================================================================
# Check it works---------
#======================================================================================
load("delay_integrals.RData")
adult_timeseries_check_mean = vector(length = length(times))
adult_timeseries_check_varying = vector(length = length(times))
adult_timeseries_check_mean[1] = adult_timeseries_trend[1]
adult_timeseries_check_varying[1] = adult_timeseries_trend[1]

for (i in seq(2, length(times))){
  adult_timeseries_check_mean[i] = adult_timeseries_check_mean[i-1] + emergence_mean[i-1]*integrate_trapezoid(integrand = emergence_integrand_mean, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_mean, lambdas) - DeathRate[i-1]*adult_timeseries_check_mean[i-1]
  adult_timeseries_check_varying[i] = adult_timeseries_check_varying[i-1] + emergence_varying[i-1]*integrate_trapezoid(integrand = emergence_integrand_varying, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_varying, egg_development_rate, larval_development_rate, pupal_development_rate) - DeathRate[i-1]*adult_timeseries_check_varying[i-1]
  print(i)
}

plot(adult_timeseries_trend)
lines(adult_timeseries_check_mean, col='red')
lines(adult_timeseries_check_varying, col='blue')


#======================================================================================
# Correct emergence timeseries for negative values, and check ---------
# 0. Do not make negative
# 1. Make negative values 0
# 2. Make negative values NA, fit a spline to the log
# 3. Take a rolling average
#======================================================================================
load("delay_integrals.RData")
emergence_mean_1 = pmax(emergence_mean, rep(0, length(emergence_mean)))
emergence_mean_2 = exp(na.spline(log(emergence_mean)))
emergence_mean_3 = emergence_mean
i=0
while (max(emergence_mean_3<0)){
  i=i+1
  print(i)
  emergence_mean_3 = rollapply(emergence_mean, i, mean, partial = TRUE)
}

emergence_varying_1 = pmax(emergence_varying, rep(0, length(emergence_varying)))
emergence_varying_2 = exp(na.spline(log(emergence_varying)))
emergence_varying_3 = emergence_varying
i=0
while (max(emergence_varying_3<0)){
  i=i+1
  print(i)
  emergence_varying_3 = rollapply(emergence_varying, i, mean, partial = TRUE)
}
adult_timeseries_check_mean_1 = vector(length = length(times))
adult_timeseries_check_mean_1[1] = adult_timeseries_trend[1]
adult_timeseries_check_mean_2 = vector(length = length(times))
adult_timeseries_check_mean_2[1] = adult_timeseries_trend[1]
adult_timeseries_check_mean_3 = vector(length = length(times))
adult_timeseries_check_mean_3[1] = adult_timeseries_trend[1]
adult_timeseries_check_varying_1 = vector(length = length(times))
adult_timeseries_check_varying_1[1] = adult_timeseries_trend[1]
adult_timeseries_check_varying_2 = vector(length = length(times))
adult_timeseries_check_varying_2[1] = adult_timeseries_trend[1]
adult_timeseries_check_varying_3 = vector(length = length(times))
adult_timeseries_check_varying_3[1] = adult_timeseries_trend[1]
for (i in seq(2, length(times))){
  adult_timeseries_check_mean_1[i] = adult_timeseries_check_mean_1[i-1] + emergence_mean_1[i-1]*integrate_trapezoid(integrand = emergence_integrand_mean, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_mean_1, lambdas) - DeathRate[i-1]*adult_timeseries_check_mean_1[i-1]
  adult_timeseries_check_mean_2[i] = adult_timeseries_check_mean_2[i-1] + emergence_mean_2[i-1]*integrate_trapezoid(integrand = emergence_integrand_mean, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_mean_2, lambdas) - DeathRate[i-1]*adult_timeseries_check_mean_2[i-1]
  adult_timeseries_check_mean_3[i] = adult_timeseries_check_mean_3[i-1] + emergence_mean_3[i-1]*integrate_trapezoid(integrand = emergence_integrand_mean, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_mean_3, lambdas) - DeathRate[i-1]*adult_timeseries_check_mean_3[i-1]
  adult_timeseries_check_varying_1[i] = adult_timeseries_check_varying_1[i-1] + emergence_varying_1[i-1]*integrate_trapezoid(integrand = emergence_integrand_varying, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_varying_1, egg_development_rate, larval_development_rate, pupal_development_rate) - DeathRate[i-1]*adult_timeseries_check_varying_1[i-1]
  adult_timeseries_check_varying_2[i] = adult_timeseries_check_varying_2[i-1] + emergence_varying_2[i-1]*integrate_trapezoid(integrand = emergence_integrand_varying, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_varying_2, egg_development_rate, larval_development_rate, pupal_development_rate) - DeathRate[i-1]*adult_timeseries_check_varying_2[i-1]
  adult_timeseries_check_varying_3[i] = adult_timeseries_check_varying_3[i-1] + emergence_varying_3[i-1]*integrate_trapezoid(integrand = emergence_integrand_varying, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_varying_3, egg_development_rate, larval_development_rate, pupal_development_rate) - DeathRate[i-1]*adult_timeseries_check_varying_3[i-1]
  print(i)
}

save(emergence_mean_1, emergence_mean_2, emergence_mean_3, emergence_varying_1, emergence_varying_2, emergence_varying_3, file = "nonneg_delay_emergence.RData")
save(adult_timeseries_check_mean, adult_timeseries_check_varying,adult_timeseries_check_mean_1, adult_timeseries_check_varying_1,adult_timeseries_check_mean_2, adult_timeseries_check_varying_2,adult_timeseries_check_mean_3, adult_timeseries_check_varying_3, file = "delayed_checks.RData")
#As expected varying_3 is the best, but a correction factor is needed
#======================================================================================
# estimate correction_factor---------
#======================================================================================
correction_factors = seq(0.991, 0.998, 0.001)
load("delay_integrals.RData")
load("nonneg_delay_emergence.RData")
adult_timeseries_check_cf_matrix = matrix(nrow = length(times), ncol=length(correction_factors))
adult_timeseries_check_cf_matrix[1,] = adult_timeseries_trend[1]

for (j in seq(length(correction_factors))){
  correction_factor = correction_factors[j]
  for (i in seq(2, length(times))){
    adult_timeseries_check_cf_matrix[i,j] = adult_timeseries_check_cf_matrix[i-1,j] + correction_factor*emergence_varying_3[i-1]*integrate_trapezoid(integrand = emergence_integrand_varying, lower = 0, upper = maxn, nstrips = number_of_strips, i-1, adult_timeseries_check_cf_matrix[,j], egg_development_rate, larval_development_rate, pupal_development_rate) - DeathRate[i-1]*adult_timeseries_check_cf_matrix[i-1,j]
  }
  print(j)
}
abs_errors=vector(length=length(correction_factors))
for (j in seq(length(correction_factors))){
  abs_errors[j] = sum(abs(adult_timeseries_trend - adult_timeseries_check_cf_matrix[,j]))
}

save(adult_timeseries_check_cf_matrix, file="cf_check.RData")
#somewhere between 0.990 and 0.999...
#best is 0.995
correction_factor = correction_factors[which(abs_errors==min(abs_errors))]
#======================================================================================
# Write to output---------
#======================================================================================
# First, the location varying parameters:
Locations$initial_adults = initial_abundances
save_columns = c("xcor","ycor","landuse","block","zone","neighborhood","code","initial_adults")
LocationsMosquitoes = Locations[save_columns]
write.csv(LocationsMosquitoes,"../output/Locations20160720.csv",row.names = F, quote = FALSE)
# save(Emergence_TS,file = sprintf("../data/Emergence%s.RData",o_f))

fileName <- file("../output/fixed_parameters_delay.txt")
open(fileName, "wt")
writeLines(paste("maxn = ", maxn), fileName)
close(fileName)


# The time-varying parameters are done in generate_parameters_climate. Save them here:
# Is the variable timeseries with method 3 the best? Check.
save(emergence_varying_3, egg_development_rate, larval_development_rate, pupal_development_rate, DeathRate, density_integral_varying, correction_factor, file="time_dependent_series_delay.RData")

