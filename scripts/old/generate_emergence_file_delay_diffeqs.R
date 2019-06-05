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
# Extend dvec 180 days in both directions, as we need a burn in period
# (period of interest 2000 - 2010-07-04 inclusive)
#======================================================================================
burn_period = 180
dvec = seq(as.Date("2000-01-01")-burn_period,as.Date("2010-12-31")+burn_period,by="1 day")
tvec = seq(as.Date("2000-01-01")-burn_period,as.Date("2010-12-31")+burn_period,by="1 day")
# tvec = seq(as.POSIXct("2000-01-01"),as.POSIXct("2010-12-31 23:00"),by="hours")
initial_date = as.Date("2000-01-01") #this is the date we normalize by, i.e. first day of interest
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
normalization_index = which(dvec == initial_date)
normalized_daily_adults = exp(daily_abundance_surf[[1]])[1,]/exp(daily_abundance_surf[[1]])[1,normalization_index]
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
# GEt pdfs etc.
#======================================================================================
load("surfaces.RData")
load("non_homogeneous_exponential_densities.RData")
times=seq(length(tvec))

alpha = 
beta = 0.571

Lambda_E = Lambda(times, egg_development_rate)
Lambda_L = Lambda(times, larval_development_rate)
Lambda_P = Lambda(times, pupal_development_rate)

f_E = dnonhomexp(times, egg_development_rate, Lambda_E)   
f_L = dnonhomexp(times, larval_development_rate, Lambda_L)   
f_P = dnonhomexp(times, pupal_development_rate, Lambda_P)   

#Estimate an nmax - just use burn period instead
# nmax=0
# acceptibility_threshold=1e-11
# for (j in times){
#   i=j
#   while ((f_E[i,j]>acceptibility_threshold | f_L[i,j]>acceptibility_threshold | f_P[i,j]>acceptibility_threshold) & i<length(times)){
#    i=i+1
#   }
#   nmax = max(nmax, i-j-1)
# }
#nmax is 72 for 1e-5, 89 for 1e-6, 107 for 1e-7, 124 for 1e-8, 141 for 1e-9, 159 for 1e-10

f_M = convolution(f_E, f_L, burn_period)
f_N = convolution(f_M, f_P, burn_period)

save(f_E, f_L, f_P, f_M, f_N, file = "non_homogeneous_exponential_densities.RData")
#======================================================================================
# GEt timeseries
# The method is not accurate for the first (or last) nmax = 150 days as need a lead time to get the integral
#======================================================================================
load("surfaces.RData")
load("non_homogeneous_exponential_densities.RData")
dNdt_with_burn_in = differentiate_timeseries(adult_timeseries_trend[1:(length(tvec)-burn_period)], tvec[1:(length(tvec)-burn_period)]) #check how this comes out - do the final bits mess it up
dNdt = dNdt_with_burn_in[(burn_period+1):(length(tvec)-burn_period)]

emergence_integral = vector(length=length(dNdt))
for (i in seq(length(dNdt))) {
  product = (adult_timeseries_trend*f_N[i+burn_period,])[1:(length(tvec)-burn_period)] #here we include the early time-points for the integral
  emergence_integral[i] = sum(product[1:(burn_period+i-1)] + product[2:(burn_period+i)])/2
}

emergence = (dNdt + (DeathRate*adult_timeseries_trend)[(burn_period+1):(length(tvec)-burn_period)])/emergence_integral

save(emergence_integral, emergence, file="delay_integrals.RData")

#======================================================================================
# Check it works---------
# The method is not accurate for the first burn_period days as need a lead time to get the integral
#======================================================================================
load("delay_integrals.RData")
adult_timeseries_check = vector(length = length(emergence)+burn_period) #needs the initial burn in period for the first timesteps
adult_timeseries_check[1:(burn_period+1)] = adult_timeseries_trend[1:(burn_period+1)]

for (i in seq(2+burn_period, length(emergence)+burn_period)){
  product = adult_timeseries_check*f_N[i-1,1:(length(tvec)-burn_period)]
  emergence_integral_check = sum(product[1:(i-1)] + product[2:i])/2
  adult_timeseries_check[i] = adult_timeseries_check[i-1] + emergence[i-burn_period-1]*emergence_integral_check - DeathRate[i-1]*adult_timeseries_check[i-1]
  print(i)
}

plot(adult_timeseries_trend[(burn_period+1):length(adult_timeseries_check)])
lines(adult_timeseries_check[(burn_period+1):length(adult_timeseries_check)], col='red')


#======================================================================================
# Correct emergence timeseries for negative values, and check ---------
# 0. Do not make negative
# 1. Make negative values 0
# 2. Make negative values NA, fit a spline to the log
# 3. Take a rolling average
#======================================================================================
load("delay_integrals.RData")
emergence_1 = pmax(emergence, rep(0, length(emergence)))
emergence_2 = exp(na.spline(log(emergence)))
emergence_3 = emergence
i=0
while (max(emergence_3<0)){
  i=i+1
  print(i)
  emergence_3 = rollapply(emergence, i, mean, partial = TRUE)
}

adult_timeseries_check_1 = vector(length = length(emergence)+burn_period)
adult_timeseries_check_1[1:(burn_period+1)] = adult_timeseries_trend[1:(burn_period+1)]
adult_timeseries_check_2 = adult_timeseries_check_1
adult_timeseries_check_3 = adult_timeseries_check_1

for (i in seq(2+burn_period, length(emergence)+burn_period)){
  product_1 = adult_timeseries_check_1*f_N[i-1,1:(length(tvec)-burn_period)]
  emergence_integral_check_1 = sum(product_1[1:(i-1)] + product_1[2:i])/2
  adult_timeseries_check_1[i] = adult_timeseries_check_1[i-1] + emergence_1[i-burn_period-1]*emergence_integral_check_1 - DeathRate[i-1]*adult_timeseries_check_1[i-1]
  product_2 = adult_timeseries_check_2*f_N[i-1,1:(length(tvec)-burn_period)]
  emergence_integral_check_2 = sum(product_2[1:(i-1)] + product_2[2:i])/2
  adult_timeseries_check_2[i] = adult_timeseries_check_2[i-1] + emergence_2[i-burn_period-1]*emergence_integral_check_2 - DeathRate[i-1]*adult_timeseries_check_2[i-1]
  product_3 = adult_timeseries_check_3*f_N[i-1,1:(length(tvec)-burn_period)]
  emergence_integral_check_3 = sum(product_3[1:(i-1)] + product_3[2:i])/2
  adult_timeseries_check_3[i] = adult_timeseries_check_3[i-1] + emergence_3[i-burn_period-1]*emergence_integral_check_3 - DeathRate[i-1]*adult_timeseries_check_3[i-1]
  print(i)
}

save(emergence_1, emergence_2, emergence_3, file = "nonneg_delay_emergence.RData")
save(adult_timeseries_check,adult_timeseries_check_1, adult_timeseries_check_2, adult_timeseries_check_3, file = "delayed_checks.RData")
#As expected varying_3 is the best, but a correction factor is needed
#======================================================================================
# estimate correction_factor---------
#======================================================================================
correction_factors = seq(0.98, 1.00, 0.001)
load("delay_integrals.RData")
load("nonneg_delay_emergence.RData")
adult_timeseries_check_cf_matrix = matrix(data=0,nrow = length(emergence)+burn_period, ncol=length(correction_factors))
adult_timeseries_check_cf_matrix[1:(burn_period+1),] = adult_timeseries_trend[1:(burn_period+1)]

for (j in seq(length(correction_factors))){
  correction_factor = correction_factors[j]
  for (i in seq(2+burn_period, length(emergence)+burn_period)){
    product = adult_timeseries_check_cf_matrix[,j]*f_N[i-1,1:(length(tvec)-burn_period)]
    emergence_integral_check = sum(product[1:(i-1)] + product[2:i])/2
    adult_timeseries_check_cf_matrix[i,j] = adult_timeseries_check_cf_matrix[i-1,j] + correction_factor*emergence_3[i-burn_period-1]*emergence_integral_check - DeathRate[i-1]*adult_timeseries_check_cf_matrix[i-1,j]
  }
  print(j)
}

abs_errors=vector(length=length(correction_factors))
for (j in seq(length(correction_factors))){
  abs_errors[j] = sum(abs(adult_timeseries_trend[(burn_period+1):(length(tvec)-burn_period)] - adult_timeseries_check_cf_matrix[(burn_period+1):(length(adult_timeseries_check_cf_matrix[,j])),j]))
}

save(adult_timeseries_check_cf_matrix, file="cf_check.RData")
#somewhere between 0.980 and 0.999...
#best is 0.989
correction_factor = correction_factors[which(abs_errors==min(abs_errors))]
plot(adult_timeseries_trend[(burn_period+1):(length(tvec)-burn_period)])
lines(adult_timeseries_check_cf_matrix[(burn_period+1):(length(adult_timeseries_check_cf_matrix[,which(abs_errors==min(abs_errors))])),which(abs_errors==min(abs_errors))], col='red')
#======================================================================================
# Write to output---------
#======================================================================================
# First, the location varying parameters:
Locations$initial_adults = initial_abundances
save_columns = c("xcor","ycor","landuse","block","zone","neighborhood","code","initial_adults")
LocationsMosquitoes = Locations[save_columns]
write.csv(LocationsMosquitoes,"../output/Locations20160720.csv",row.names = F, quote = FALSE)
# save(Emergence_TS,file = sprintf("../data/Emergence%s.RData",o_f))

densityFileName <- "density_values.csv"
integration_period = burn_period
f_N_output = matrix(0,ncol=integration_period, nrow=length(emergence))
f_E_output =f_N_output
f_L_output =f_N_output
f_P_output =f_N_output
f_M_output =f_N_output
for (i in seq(dim(f_N_output)[1])){
  f_N_output[i,] = f_N[i+burn_period,(i+1+burn_period-integration_period):(i+burn_period)]
  f_E_output[i,] = f_E[i+burn_period,(i+1+burn_period-integration_period):(i+burn_period)]
  f_L_output[i,] = f_L[i+burn_period,(i+1+burn_period-integration_period):(i+burn_period)]
  f_P_output[i,] = f_P[i+burn_period,(i+1+burn_period-integration_period):(i+burn_period)]
  f_M_output[i,] = f_M[i+burn_period,(i+1+burn_period-integration_period):(i+burn_period)]
}
write.table(f_N_output,paste0("../output/", densityFileName),row.names = F, col.names = F,quote = FALSE, sep=',')
write.table(f_E_output,"../output/egg_density.csv",row.names = F, col.names = F,quote = FALSE, sep=',')
write.table(f_L_output,"../output/larval_density.csv",row.names = F, col.names = F,quote = FALSE, sep=',')
write.table(f_P_output,"../output/pupal_density.csv",row.names = F, col.names = F,quote = FALSE, sep=',')
write.table(f_M_output,"../output/middle_density.csv",row.names = F, col.names = F,quote = FALSE, sep=',')

popnFileName <- "burnin_population.csv"
write.table(adult_timeseries_trend[(burn_period-integration_period+1):burn_period],paste0("../output/",popnFileName),row.names = F, col.names = F,quote = FALSE, sep=',')

fileName <- file("../output/fixed_parameters_delay.txt")
open(fileName, "wt")
writeLines(paste("burnin_population = ", popnFileName), fileName)
writeLines(paste("density_values = ", densityFileName), fileName)
close(fileName)


# The time-varying parameters are done in generate_parameters_climate. Save them here:
# Is the variable timeseries with method 3 the best? Check.
DeathRate = DeathRate[(1+burn_period):(length(tvec)-burn_period)]
save(emergence_3, DeathRate, correction_factor, file="time_dependent_series_delay.RData")

