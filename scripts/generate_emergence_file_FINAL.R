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
               zoo,
               doParallel,
               pspline,
               viridis,
               deSolve,
               scam,
               mgcv,
               here)

# Set your working directory to the scripts folder
setwd(here())
setwd("/home/sean/Documents/zika_project/mosquito_dynamics/scripts/")

source('./generate_emergence_functions_FINAL.R')
source('response_curves.R')


#======================================================================================
# Load data -----------------
#======================================================================================
dvec = seq(as.Date("2000-01-01"),as.Date("2010-12-31"),by="1 day")
tvec = seq(as.Date("2000-01-01"),as.Date("2010-12-31"),by="1 day")
# tvec = seq(as.POSIXct("2000-01-01"),as.POSIXct("2010-12-31 23:00"),by="hours")
initial_date = as.Date("2000-01-01")
# Locations = read.csv("../data/locations_area_20140807.csv")
Locations = read.csv("../data/synthetic_locations.csv")
Locations_original = Locations
# House_area = mean(Locations_original$area[which(Locations_original$landuse=="HOUSE")])


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

# Locations$ProdRatioArea = ProdRatio * Locations$area /  House_area
# Locations$ProdRatioArea[which(Locations$landuse == "HOUSE")] = 1
# Locations$ProdRatioArea[which(Locations$ProdRatioArea > 2)] = 2

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

## Plot the temperature time series (mean and max)
pdf("../figures/temperature_time_series.pdf")
par(mfrow=c(2,3),mar=c(4.1,4.1,1.1,1.1))
plot(seq(as.Date("2000-01-01"),length.out=length(temperature),by="1 day"),temperature,
     xlab="Year",ylab=expression("Daily mean air temperature ("*~degree*C*")"), las=1, xaxs="i",yaxs="i", bty="n",type='l')
plot(seq(as.Date("2000-01-01"),length.out=length(temperature.max),by="1 day"),temperature.max,
     xlab="Year",ylab=expression("Daily maximum air temperature ("*~degree*C*")"), las=1, xaxs="i",yaxs="i", bty="n",type='l')
plot(seq(as.Date("2000-01-01"),length.out=length(temperature.min),by="1 day"),temperature.min,
     xlab="Year",ylab=expression("Daily minimum air temperature ("*~degree*C*")"), las=1, xaxs="i",yaxs="i", bty="n",type='l')
plot(seq(as.Date("2000-01-01"),length.out=length(water.temp.mean),by="1 day"),water.temp.mean,
     xlab="Year",ylab=expression("Daily mean water temperature ("*~degree*C*")"), las=1, xaxs="i",yaxs="i", bty="n",type='l')
plot(seq(as.Date("2000-01-01"),length.out=length(water.temp.max),by="1 day"),water.temp.max,
     xlab="Year",ylab=expression("Daily maximum water temperature ("*~degree*C*")"), las=1, xaxs="i",yaxs="i", bty="n",type='l')
plot(seq(as.Date("2000-01-01"),length.out=length(water.temp.min),by="1 day"),water.temp.min,
     xlab="Year",ylab=expression("Daily minimum water temperature ("*~degree*C*")"), las=1, xaxs="i",yaxs="i", bty="n",type='l')
dev.off()
#======================================================================================
# obtain parameters governing death and development ------------
#======================================================================================
development_params = read.csv('../data/development_params.csv', sep = '\t', row.names=1)
mortality_thresholds = read.csv('../data/mortality_thresholds.csv', sep = '\t', row.names=1)

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
# spray_2004_mid = dvec[1821]
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
#log_spray_abundances = Estimate.w.offset(XYall,spray_2004_mid,as.numeric(o_f))
#spray_abundances = proto_pac*exp(log_spray_abundances[[1]])[,1]
#Some of the new points are outside the old polygon - choose the nearest point and assign that initial abundance.
non_na_coordinates = which(!is.na(initial_abundances))
non_na_initial_abundances = initial_abundances[non_na_coordinates]
#non_na_spray_abundances = spray_abundances[non_na_coordinates]
for (i in which(is.na(initial_abundances))){
  distances = sqrt((Locations$xcor[non_na_coordinates] - Locations$xcor[i])^2 + (Locations$ycor[non_na_coordinates] - Locations$ycor[i])^2)
  initial_abundances[i] = non_na_initial_abundances[which.min(distances)]
#  spray_abundances[i] = non_na_spray_abundances[which.min(distances)]
}

#sprayed:
#daily_abundance_surf.sprayed = Estimate.w.offset.sprayed(XY,dvec,as.numeric(o_f))
#normalized_daily_adults = exp(daily_abundance_surf.sprayed[[1]])[1,]/exp(daily_abundance_surf.sprayed[[1]])[1,1]
#normalized_adults = rep(NA, length=length(tvec))
#normalized_adults[seq(floor(length(tvec)/length(dvec)/2+1), length(tvec), length(tvec)/length(dvec))] <- normalized_daily_adults
#adult_timeseries_trend.sprayed = na.spline(normalized_adults)
#log_initial_abundances = Estimate.w.offset.sprayed(XYall,initial_date,as.numeric(o_f))
#initial_abundances.sprayed = proto_pac*exp(log_initial_abundances[[1]])[,1]
#log_spray_abundances = Estimate.w.offset.sprayed(XYall,spray_2004_mid,as.numeric(o_f))
#spray_abundances.sprayed = proto_pac*exp(log_spray_abundances[[1]])[,1]
#Some of the new points are outside the old polygon - choose the nearest point and assign that initial abundance.
#non_na_coordinates = which(!is.na(initial_abundances.sprayed))
#non_na_initial_abundances = initial_abundances.sprayed[non_na_coordinates]
#non_na_spray_abundances = spray_abundances.sprayed[non_na_coordinates]
#for (i in which(is.na(initial_abundances.sprayed))){
#  distances = sqrt((Locations$xcor[non_na_coordinates] - Locations$xcor[i])^2 + (Locations$ycor[non_na_coordinates] - Locations$ycor[i])^2)
#  initial_abundances.sprayed[i] = non_na_initial_abundances[which.min(distances)]
#  spray_abundances.sprayed[i] = non_na_spray_abundances[which.min(distances)]
#}
#save(adult_timeseries_trend, initial_abundances, daily_abundance_surf,adult_timeseries_trend.sprayed, initial_abundances.sprayed, daily_abundance_surf.sprayed, file="early_surfaces.RData")
save(adult_timeseries_trend, initial_abundances, daily_abundance_surf, file="early_surfaces.RData")
#To get the trend for location i, do initial_abundances[i]*trend_abundance_surf

write.csv(adult_timeseries_trend*sum(initial_abundances), "../output/mosquito_ts.csv", row.names = F, quote = F)

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

cols <- c(red="#AE1A19", black="#252525", grey="#B1B1A6",
          brown="#604B4A", lgrey="#D0D0D0", white="#FFFFFF")

larval_capacity = 0.2 #set to Inf for no density dependence
larval_power = 2
eggs_per_gon_cycle = 63

##get the simple one-step model first.
emergence.ts <- differentiate_timeseries(adult_timeseries_trend, tvec) + DeathRate*adult_timeseries_trend
emergence.ts.nn1 = pmax(emergence.ts, rep(0, length(emergence.ts)))
emergence.ts.nn2 = exp(na.spline(log(emergence.ts)))
emergence.ts.nn3 = emergence.ts
i=0
while (max(emergence.ts.nn3<0)){
  i=i+1
  print(i)
  emergence.ts.nn3 = rollapply(emergence.ts, i, mean, partial = TRUE)
}
emergence.pop = matrix(data=0, nrow=length(adult_timeseries_trend), ncol=3)
parms1 <- list(emergence.ts.nn1, DeathRate)
parms2 <- list(emergence.ts.nn2, DeathRate)
parms3 <- list(emergence.ts.nn3, DeathRate)
emergence.pop[,1] = ode(y = c(N=adult_timeseries_trend[1]), times = seq(1,length(tvec),1), func = emergence.ode, parms = parms1)[,"N"]
emergence.pop[,2] = ode(y = c(N=adult_timeseries_trend[1]), times = seq(1,length(tvec),1), func = emergence.ode, parms = parms2)[,"N"]
emergence.pop[,3] = ode(y = c(N=adult_timeseries_trend[1]), times = seq(1,length(tvec),1), func = emergence.ode, parms = parms3)[,"N"]

abs_errors = vector(length=3)
for (i in 1:3){
  abs_errors[i] = sum(abs(adult_timeseries_trend-emergence.pop[,i])/adult_timeseries_trend)
}

minimize_function = function(correction_factor, DeathRate, emergence.ts.nn3){
    params = list(correction_factor*emergence.ts.nn3, DeathRate)
    population_timeseries = ode(y = c(N=adult_timeseries_trend[1]), times = seq(1,length(tvec),1), func = emergence.ode, parms = params)
  return(sum(abs(adult_timeseries_trend-population_timeseries[,"N"])/adult_timeseries_trend))
}

corrfactor.optimized.emergence = optimize(minimize_function, DeathRate,
                                          emergence.ts.nn3, interval=c(0.8,1.25))

corrfactor.emergence = corrfactor.optimized.emergence$minimum
corrfactor.emergence <- 0.9880937

##now the full model, withe early stages
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
pupal_test_pop = ode(y = c(N=1), times = seq(1,length(tvec),1), func = pupal_test_ode, parms = list(pupal_development_rate, pupal_timeseries_nn6, DeathRate), method='radau')[,"N"]
save(pupal_timeseries, pupal_timeseries_nn2, pupal_timeseries_nn4, pupal_timeseries_nn6, file="pupal_timeseries.RData")


egg_timeseries = get_egg_timeseries_ode_2(adult_timeseries_trend,egg_mortality, egg_development_rate,gonotrophic_cycle_rate,eggs_per_gon_cycle, 12, tvec, method='radau')
print(min(egg_timeseries))
save(egg_timeseries, file="egg_timeseries.RData")


#larval_timeseries_nn2 = get_larval_timeseries_ode(egg_surf.in = egg_timeseries, pupal_surf.in = pupal_timeseries_nn2, egg_development_rate.in = egg_development_rate, larval_development_rate.in = larval_development_rate, pupal_development_rate.in = pupal_development_rate, larval_mortality.in = larval_mortality, pupal_mortality.in = pupal_mortality, larval_capacity.in = larval_capacity, larval_power.in = larval_power, initial_larvae = 1, times.in = tvec)
larval_timeseries_nn4 = get_larval_timeseries_ode(egg_surf.in = egg_timeseries, pupal_surf.in = pupal_timeseries_nn4, egg_development_rate.in = egg_development_rate, larval_development_rate.in = larval_development_rate, pupal_development_rate.in = pupal_development_rate, larval_mortality.in = larval_mortality, pupal_mortality.in = pupal_mortality, larval_capacity.in = larval_capacity, larval_power.in = larval_power, initial_larvae = 0.66, times.in = tvec, method='radau')
larval_timeseries_nn6 = get_larval_timeseries_ode(egg_surf.in = egg_timeseries, pupal_surf.in = pupal_timeseries_nn6, egg_development_rate.in = egg_development_rate, larval_development_rate.in = larval_development_rate, pupal_development_rate.in = pupal_development_rate, larval_mortality.in = larval_mortality, pupal_mortality.in = pupal_mortality, larval_capacity.in = larval_capacity, larval_power.in = larval_power, initial_larvae = 0.66, times.in = tvec, method='radau')
print(min(larval_timeseries_nn4, larval_timeseries_nn6))
save(larval_timeseries_nn4, larval_timeseries_nn6, file="larval_timeseries.RData")
##New method does not work without using non-negative and non-zero timeseries - i.e. only methods 6 and 4 work.


extra_mortality_nn4 = (egg_development_rate*egg_timeseries - differentiate_timeseries(larval_timeseries_nn6, tvec))/larval_timeseries_nn6 - larval_mortality - larval_development_rate - larval_timeseries_nn4^larval_power/larval_capacity
extra_mortality_nn4 = exp(na.spline(log(extra_mortality_nn4)))
extra_mortality_nn6 = (egg_development_rate*egg_timeseries - differentiate_timeseries(larval_timeseries_nn6, tvec))/larval_timeseries_nn6 - larval_mortality - larval_development_rate - larval_timeseries_nn6^larval_power/larval_capacity
print(min(extra_mortality_nn6))
extra_mortality_nn6_initial = extra_mortality_nn6 
i=0
while (max(extra_mortality_nn6 <0)){
  i=i+1
  print(i)
  extra_mortality_nn6 = rollapply(extra_mortality_nn6_initial, i, mean, partial = TRUE)
}

##plot mortalities
pupal.mort.summary <- data.frame("envt"=pupal_mortality, "extr"=extra_mortality_nn6)
tot.pupal.mort <- apply(pupal.mort.summary, 1, sum)
larval.mort.summary <- data.frame("envt"=pupal_mortality, "extr"=extra_mortality_nn6,
                                  "dens"=larval_timeseries_nn6^larval_power/larval_capacity)
tot.larval.mort <- apply(larval.mort.summary, 1, sum)
cols <- plasma(3)

tiff(paste0("../figures/pupal_mortality_source.tif"),
     res=600, width=4152, height=4152, compression="lzw")
plot(x=dvec, y=rep(-1, length(dvec)),
     col=cols[1], type='l', ylim=c(0,1),
     xlab="Year", ylab="Proportion")
polygon(x=c(dvec, rev(dvec)),
        y=c(rep(0, length(dvec)), rev(pupal.mort.summary$extr/tot.pupal.mort)),
        border=NA, col=cols[1])
polygon(x=c(dvec, rev(dvec)),
        y=c(pupal.mort.summary$extr/tot.pupal.mort, rev((pupal.mort.summary$envt+pupal.mort.summary$extr)/tot.pupal.mort)),
        border=NA, col=cols[2])
legend(x=as.Date("2005-09-01"), y=0.58, legend=c("Temperature driven", "Additional"),
       col=cols[c(2, 1)], fill=cols[c(2, 1)])
dev.off()

tiff(paste0("../figures/larval_mortality_source.tif"),
     res=600, width=4152, height=4152, compression="lzw")
plot(x=dvec, y=rep(-1, length(dvec)),
     col=cols[1], type='l', ylim=c(0,1),
     xlab="Year", ylab="Proportion")
polygon(x=c(dvec, rev(dvec)),
        y=c(rep(0, length(dvec)), rev(larval.mort.summary$extr/tot.larval.mort)),
        border=NA, col=cols[1])
polygon(x=c(dvec, rev(dvec)),
        y=c(larval.mort.summary$extr/tot.larval.mort, rev((larval.mort.summary$extr+larval.mort.summary$envt)/tot.larval.mort)),
        border=NA, col=cols[2])
polygon(x=c(dvec, rev(dvec)),
        y=c((larval.mort.summary$extr+larval.mort.summary$envt)/tot.larval.mort, rep(1, length(dvec))),
        border=NA, col=cols[3])
legend(x=as.Date("2005-09-01"), y=0.99, legend=c( "Density dependent","Temperature driven", "Additional"),
       col=cols[3:1], fill=cols[3:1])
dev.off()

## check_larvae = (differentiate_timeseries(pupal_timeseries_nn6, tvec) + (pupal_mortality+extra_mortality_nn6+pupal_development_rate)*pupal_timeseries_nn6)/larval_development_rate
# summary(larval_timeseries_nn6/check_larvae)
# check_extra_mortality = (larval_development_rate*larval_timeseries_nn6 - differentiate_timeseries(pupal_timeseries_nn6, tvec))/pupal_timeseries_nn6 - larval_mortality - pupal_mortality
# summary(extra_mortality_nn6/check_extra_mortality)

# Methods to compare========================================
# 4. Make negative values NA, fit a spline to the log, thoughout
# 6. Take a rolling average throughout
#===========================================================
population_timeseries = matrix(data=0, nrow=length(adult_timeseries_trend), ncol=2)
initial_conditions_4 = c(E =egg_timeseries[1], L = larval_timeseries_nn4[1], P = pupal_timeseries_nn4[1], N = adult_timeseries_trend[1])
initial_conditions_6 = c(E =egg_timeseries[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])
params_4 = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn4, DeathRate, larval_capacity, larval_power, eggs_per_gon_cycle)
params_6 = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, eggs_per_gon_cycle)
#params_check = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, check_extra_mortality, DeathRate, larval_capacity, larval_power, eggs_per_gon_cycle)
population_timeseries[,1] = ode(y = initial_conditions_4, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_4, method='radau')[,"N"]
population_timeseries[,2] = ode(y = initial_conditions_6, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_6, method='radau')[,"N"]
#population_timeseries[,3] = ode(y = initial_conditions_6, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_check)[,"N"]

abs_errors = vector(length=2)
for (i in 1:2){
  abs_errors[i] = sum(abs(adult_timeseries_trend-population_timeseries[,i])/adult_timeseries_trend)
}

###correction factor is not necessaryily needed - but see what fits best anyways
initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])

minimize_function = function(correction_factor, DeathRate, pupal_mortality, larval_mortality, egg_mortality, extra_mortality_nn6, adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity, larval_power, eggs_per_gon_cycle){
  params = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, correction_factor*eggs_per_gon_cycle)
  population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params)
  return(sum(abs(adult_timeseries_trend-population_timeseries[,"N"])/adult_timeseries_trend))
}

correction_factor_optimized = optimize(minimize_function, DeathRate, pupal_mortality, larval_mortality, egg_mortality, extra_mortality_nn6, adult_development_rate, pupal_development_rate, larval_development_rate, egg_development_rate,  larval_capacity, larval_power, eggs_per_gon_cycle, interval=c(0.8,1.25))

correction_factor = correction_factor_optimized$minimum

correction_factor <- 1.056433
params = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, correction_factor*eggs_per_gon_cycle)
#population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params)


##import abm output
abm.out <- read.csv("../scripts/sim_00020_immature.csv")["Adults"]
abm.out.irs <- read.csv("../scripts/sim_00010_immature.csv")["Adults"]
abm.out.ulv <- read.csv("../scripts/sim_00000_immature.csv")["Adults"]
cols <- plasma(4, alpha=c(0.2,rep(1,2)))

tiff("../figures/abundance_vs_time.tif",
     res=600, width=4152, height=4152/2, compression="lzw")
par(mar=c(5.1,5.1,4.1,2.1))
plot(dvec, adult_timeseries_trend*sum(initial_abundances), ylim = c(0,1.8*sum(initial_abundances)),
     col=cols[1], type="l",lwd=2,xaxs="i",yaxs="i",bty="n",las=1,xlim=c(as.Date("2000-01-01"),as.Date("2011-01-01")),
     xlab="Year", ylab="",xaxt="n")
axis(1,label=seq(2000,2011,1),at=seq(as.Date("2000-01-01"),as.Date("2011-01-01"),by="1 year"))
mtext("Abundance",2,4)
lines(dvec,population_timeseries[,2]*sum(initial_abundances), col=cols[2], lwd=2)
lines(dvec,abm.out[,1], col=cols[3],  lwd=2)
legend("topright", legend=c("GAM prediction","ODE calibration", "ABM output"),
       col=cols[c(1,2,3)], lty = "solid",
       lwd=c(2),bty="n")
dev.off()

print(cor(adult_timeseries_trend*sum(initial_abundances),population_timeseries[,2]*sum(initial_abundances)))
print(cor(adult_timeseries_trend*sum(initial_abundances),abm.out[,1]))
print(cor(population_timeseries[,2]*sum(initial_abundances),abm.out[,1]))

##RMSD
rmsd.ode <- sum(initial_abundances)*(sum((adult_timeseries_trend-population_timeseries[,2])^2)/length(adult_timeseries_trend))^0.5
rmsd.abm <- (sum((sum(initial_abundances)*adult_timeseries_trend-abm.out[,1])^2)/length(adult_timeseries_trend))^0.5
rmsd.odeabm <- (sum((sum(initial_abundances)*population_timeseries[,2]-abm.out[,1])^2)/length(adult_timeseries_trend))^0.5

cols.spray <- plasma(9)
sday <- 366; eday <- 365*4
#ABM spraying!
tiff("../figures/abundance_following_spraying.tif",
     res=600, width=4152, height=4152/2, compression="lzw")
par(mar=c(5.1,5.1,4.1,2.1))
plot(dvec[sday:eday], abm.out.ulv[sday:eday,1], ylim = c(0,1.8*sum(initial_abundances)),
     col=cols.spray[1], las=1, bty="n",xaxs="i",yaxs="i",
     xlab="Year", ylab="",type='l', lwd=2,
     xlim=c(as.Date("2001-01-01"),as.Date("2004-01-01")))
lines(dvec[sday:eday],abm.out.irs[sday:eday,1], col=cols.spray[4], lwd=2)
lines(dvec[sday:eday],abm.out[sday:eday,1], col=cols.spray[7], lwd=2)
mtext("Abundance",2,4)
legend("topright", legend=c("TIRS","ULV", "None"), col=cols.spray[c(1,4,7)],
       lty = c("solid"),lwd=c(2,2,2),bty="n")
dev.off()

##Rebound times - remember the variables are mislabeled (correct this)
(head(which(((abm.out[,1]-abm.out.irs[,1])/abm.out[,1] > 0)
           & ((abm.out[,1]-abm.out.irs[,1])/abm.out[,1] < 0.01))) - 731)*12/365
(head(which(((abm.out[,1]-abm.out.irs[,1])/abm.out[,1] > 0)
           & ((abm.out[,1]-abm.out.irs[,1])/abm.out[,1] < 0.1))) - 731)*12/365
(head(which(((abm.out[,1]-abm.out.ulv[,1])/abm.out[,1] > 0)
           & ((abm.out[,1]-abm.out.ulv[,1])/abm.out[,1] < 0.01))) - 731)*12/365
(head(which(((abm.out[,1]-abm.out.ulv[,1])/abm.out[,1] > 0)
           & ((abm.out[,1]-abm.out.ulv[,1])/abm.out[,1] < 0.1))) - 731)*12/365

tiff("../figures/residuals.tif",
     res=600, width=4152, height=4152, compression="lzw")
plot(dvec, (population_timeseries[,2]-adult_timeseries_trend)*sum(initial_abundances),
     col=cols[2],
     xlab="Year", ylab="Residual", type='l', lwd=3)
lines(dvec, abm.out[,1]-adult_timeseries_trend*sum(initial_abundances),
      col=cols[3], lwd=3, lty="dotted")
dev.off()

##==================================================
# comparing with and without spraying -------------
#==================================================
spraying_thoroughness_vec = c(0.1)
spray_today = vector(length = length(DeathRate))
start_day = 250; end_day= start_day+21
spray_today[start_day:end_day] = TRUE 

initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries_nn6[1], P = pupal_timeseries_nn6[1], N = adult_timeseries_trend[1])
params_no_spraying = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, correction_factor*eggs_per_gon_cycle)
population_timeseries_no_spraying = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_no_spraying, method = 'radau')[,"N"]

pdf("../figures/spraying_fig_ode.pdf")
plot(seq(1,length(tvec),1)[200:600], population_timeseries_no_spraying[200:600],
     col=cols['red'], lwd=1.5, ylim=c(0,1.7), type='l', xlab="Day", ylab="Normalized abundance")
for (spraying_thoroughness in spraying_thoroughness_vec) {
    params_spraying = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate+spraying_thoroughness*spray_today, larval_capacity, larval_power, correction_factor*eggs_per_gon_cycle)
    population_timeseries_spraying = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_spraying, method = 'radau')[,"N"]
    lines(seq(1,length(tvec),1)[200:600], population_timeseries_spraying[200:600], col=cols['black'], lwd=1.5)
}
abline(v=start_day, col=cols['grey'])
abline(v=end_day, col=cols['grey'])
legend(x=412, y=0.3, legend=c("Abundance, no spraying","Abundance, spraying","Start/end spraying"),
       col=cols[c("red", "black", "grey")], lty=rep("solid", 3))
dev.off()

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
params = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity, larval_power, eggs_per_gon_cycle)
population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params, method = 'radau')[,"N"]

##tiff("../figures/rebound_from_zero.tif",
     res=600, width=4152, height=4152, compression="lzw")
plot(seq(1,length(tvec),1)[1:365], population_timeseries[1:365], type='l', col=cols["red"])
points(adult_timeseries_trend, col=cols['black'])
##dev.off()
#==================================================
# does the density dependence maintain spatial heterogeneity? 
#=================================================
#Look at the extremes, start with 0 adults
maxN0 = max(initial_abundances)
minN0 = min(initial_abundances)

initial_conditions = c(E =1e-30, L = 0, P = 0, N = 0)
params_max = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity*maxN0^larval_power, larval_power, eggs_per_gon_cycle)
params_min = list(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, extra_mortality_nn6, DeathRate, larval_capacity*minN0^larval_power, larval_power, eggs_per_gon_cycle)
population_timeseries_max = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_max, method = 'radau')[,"N"]
population_timeseries_min = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode_2, parms = params_min, method = 'radau')[,"N"]

tiff("../figures/larval_capacity_equilibrium.tif",
     res=600, width=4152, height=4152, compression="lzw")
plot(tvec[1:1250], population_timeseries_max[1:1250], type='l', col=cols[2], ylim = c(0,35),
     xlab="Year", ylab="Abundance", lwd=2)
lines(tvec[1:1250], adult_timeseries_trend[1:1250]*maxN0, col=cols[1], lty='solid', lwd=2)
lines(tvec[1:1250], adult_timeseries_trend[1:1250]*minN0, col=cols[1], lty='dotted', lwd=2)
lines(tvec[1:1250], population_timeseries_min[1:1250], col=cols[2], lty="dotted", lwd=2)
legend(x=as.Date("2002-01-01"), y=35, legend=c("largest location, ODE","largest location, GAM",
                                               "smallest location, ODE","smallest location, GAM"),
       col=cols[c(2,2,1,1)], lty = c("solid", "dashed", "solid", "dashed"), lwd=rep(2,4))
dev.off()


#======================================================================================
# Write to output---------
#======================================================================================
# First, the location varying parameters:
Locations$initial_adults = initial_abundances
# save_columns = c("xcor","ycor","landuse","block","zone","neighborhood","code","initial_adults", "moh_zone")
save_columns = c("xcor","ycor","landuse","initial_adults", "moh_zone", "locID")
LocationsMosquitoes = Locations[save_columns]
write.csv(LocationsMosquitoes,"../output/Locations20190109.csv",row.names = F, quote = FALSE)

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
eggs_per_gon_cycle  <- eggs_per_gon_cycle*correction_factor
emergence.ts <- emergence.ts.nn3*corrfactor.emergence

save(egg_development_rate, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, egg_mortality, larval_mortality, pupal_mortality, DeathRate, extra_mortality_nn6, eggs_per_gon_cycle, emergence.ts, file="time_dependent_series.RData")

