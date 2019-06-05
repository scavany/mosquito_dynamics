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

source('./generate_emergence_functions.R')
source('response_curves.R')


#======================================================================================
# Load data -----------------
#======================================================================================
tvec = seq(as.Date("2000-01-01"),as.Date("2010-12-31"),by="1 day")
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

# 100 Localities are used to calculate the relationship between emergence and space
SmallLocations = Locations[10000,]  


#======================================================================================
# Read temperature specific for Iquitos ------
#======================================================================================
Iquitos.climate = read.csv('../data/Iquitos_Climate_Bobby.csv')
startIndex = which(Iquitos.climate$date == as.character(tvec[1])) 
temperature = Iquitos.climate$temperature_mean[startIndex:(startIndex + length(tvec)-1) ]
temp.mean = mean(temperature,na.rm = TRUE); temperature[is.na(temperature)] = temp.mean


#======================================================================================
# calculate mortality ------------
#======================================================================================
DeathRate = get_mortality_vector(tvec, temperature)

#======================================================================================
# Calculate Surface Space + Time for a small dataset --------
#======================================================================================
offsets = c("1.1","1.35","1.6","1.7","1.9", "2")
# small_abundance_surf = list()
o_f = offsets[6]
eval(parse(text = sprintf("load(\"EstimateSurface%s%s\")",o_f,".RData")))
XY = cbind(SmallLocations$xcor,SmallLocations$ycor)
small_abundance_surf = Estimate.w.offset(XY,tvec,as.numeric(o_f)) # Log -- Est: [[1]] Std: [[2]]
Emergence.data = get_emergence_timeseries(1,exp(small_abundance_surf[[1]]), DeathRate, tvec)
Emergence_TS = Emergence.data$Emergence_Timeseries /  Emergence.data$Emergence_Timeseries[1]

save(Emergence_TS,file = sprintf("../data/Emergence%s.RData",o_f))

#======================================================================================
# Calc. the relationship of initial  abundance and emergence rates in time -----------
#======================================================================================
emerge_correction_factor = array(0,length(SmallLocations[,1]))
Emergence.data = get_emergence_timeseries(1,exp(small_abundance_surf[[1]]), DeathRate, tvec)

# The ODE uses Emergence_TS and deathrate as global variables
# As the temporal trend is constant in time, the timeseries are proportional between locations
# As we want to know the emergence trend given the first time point, we normalise by timepoint 1
# this normalized trend is identical across locations, and whether proto_pac is used
Emergence_TS = Emergence.data$Emergence_Timeseries /  Emergence.data$Emergence_Timeseries[1] 
params.init = c(cor_factor = 1)
times = seq(1,length(DeathRate),1)
for(i in 1:length(SmallLocations[,1])){
  #Emergence.data should have the same trend for all locations
  Emergence.data =
    get_emergence_timeseries(i, exp(small_abundance_surf[[1]]) * proto_pac, DeathRate, tvec)
  # why not make Emergence.data the global variable from which epsilon is estimated?
  yinit = c(N = Emergence.data$N[1])
  out.N1 = ode(yinit,times,pop.model.ode,params.init)
  # Calculate the mismatch between the mosquito abundance data and the ODE output
  emerge_correction_factor[i] = 
    mean(Emergence.data$N[2:length(Emergence.data$N)-1] / out.N1[2:length(Emergence.data$N)-1,2])
}



# the relationship is a straight line with intercept in 0, so I'm calculating the slope and that's it !
m_emerge = (max(emerge_correction_factor) - min(emerge_correction_factor)) / 
  (max(exp(small_abundance_surf[[1]][,1]) * proto_pac) - min(exp(small_abundance_surf[[1]][,1])) * proto_pac) 


#======================================================================================
# Calculate Total Surface Space in time 0 -----------
#======================================================================================
XY = cbind(Locations$xcor,Locations$ycor)
abundance_surf_init = Estimate.w.offset(XY,tvec,as.numeric(o_f)) # Log -- Est: [[1]] Std: [[2]]

#======================================================================================
# Write to output---------
#======================================================================================
Locations$emerge = m_emerge * exp(abundance_surf_init[[1]]) * proto_pac
LocationsMosquitoes = Locations[c(1:7,11)]
write.csv(LocationsMosquitoes,"../output/Locations20160720.csv",row.names = F, quote = FALSE)

# save(Emergence_TS,file = sprintf("../data/Emergence%s.RData",o_f))
