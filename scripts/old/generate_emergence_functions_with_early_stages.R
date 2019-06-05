#=============================================================================
# Author: Guido Espana
# 25/10/2016 
# Functions for Emergence Rates Generation
#
#=============================================================================
get_mortality_vector = function(time.in, temperature.in, fieldcorxn = 0.08969818){
  death.temp = array(data = 0, dim = c(length(time.in)))
  registerDoParallel(cores=4)
  death.temp = foreach(t = 1:length(time.in), .combine = 'c') %dopar% {
    death.temp[t] = mortalityRT(temperature.in[t],fieldcorxn) # Brady etal,2013 
  }
  return(death.temp)
}

# Calculating the Population Dynamics of one location and time 
get_abundance_in_time <- function(i,locTemp){
  N_temp <-ExpansionFactor*exp(predict(
    Output2,data.frame(Multiplier=Multiplier,YearDay=SimYears+SimDays/365,Day=SimDays,           # DeathRate*
                       Xcor=rep(locTemp$xcor[i],length(SimDays)),
                       Ycor=rep(locTemp$ycor[i],length(SimDays)),
                       Effort=rep(0,length(SimDays)))))
  N_temp = N_temp * locTemp$ProdRatioArea[i]
  return(N_temp)
}

#Original method
get_emergence_timeseries = function(i,surf.in, death_rate.in, times.in){
  registerDoParallel(cores=4)
  Emergence_Timeseries = array(data =0,dim = c(length(times.in)))
  N = surf.in[i,]
  N[is.na(N)] = 0
  DeltaN = diff(N)/as.numeric(diff(times.in))
  Emergence_Timeseries = foreach(j = 2:length(N), .combine='c') %dopar% { # for each point in time, from day 2 because we can't derive delta for day 1
    Death = death_rate.in[j]*N[j]
    Emergence_Timeseries[j]  = DeltaN[j-1] + Death     # Using backward Euler
  }
  return (list(N = N, Emergence_Timeseries = Emergence_Timeseries))
}



#original method
pop.model.ode <- function (t, x, params) {
  N <- x[1];
  cor_factor <- params[1]
  mu <- DeathRate[t+1]; epsilon <- EmergenceTS[t] ;
  dN <-  epsilon * cor_factor - mu * N
  return(list(dN))
}



#=============================================================================
# Author: Sean M. Cavany
# 14/08/2018 
# Functions for Mosquito model
#
#=============================================================================

# #My adapted method
# get_emergence_timeseries = function(i,surf.in, death_rate.in, times.in){
#   Emergence_Timeseries = array(data =0,dim = c(length(times.in)))
#   N = surf.in[i,]
#   #Is it okay to just set the population to zero for NaNs?
#   N[is.na(N)] = 0
#   DeltaN = diff(N)/as.numeric(diff(times.in))
#   Emergence_Timeseries[1]=DeltaN[1]+death_rate.in[1]*N[1]
#   Emergence_Timeseries[length(N)]=DeltaN[length(N)-1]+death_rate.in[length(N)]*N[length(N)]
#   for(j in 2:(length(N)-1)) { # for each point in time, from day 2 because we can't derive delta for day 1
#     Death = death_rate.in[j]*N[j]
#     Emergence_Timeseries[j]  = (DeltaN[j] + DeltaN[j-1])/2 + Death     # Using backward Euler
#   }
#   return (list(N = N, Emergence_Timeseries = Emergence_Timeseries))
# }

# #spline method
# get_emergence_timeseries = function(i,surf.in, death_rate.in, times.in){
#   N = surf.in[i,]
#   #Is it okay to just set the population to zero for NaNs?
#   times.in.spline = times.in
#   times.in.spline[is.na(N)] = NA
#   DeltaN = predict(sm.spline(times.in.spline, N), times.in, 1)[,1]
#   Emergence_Timeseries = DeltaN + N*death_rate.in
#   return (list(N = N, Emergence_Timeseries = Emergence_Timeseries))
# }

# #either of new methods
# pop.model.ode <- function (t, x, params) {
#   N <- x[1]; 
#   cor_factor <- params[1]
#   mu <- DeathRate[t]; 
#   epsilon <- Emergence.data$Emergence_Timeseries[t] ; 
#   dN <-  epsilon*cor_factor - mu * N
#   return(list(dN))
# }

calculate_development_rate_timeseries = function(development_params.in, temperature.in){
  rho25 = development_params.in[1]
  DeltaHA = development_params.in[2]
  DeltaHH = development_params.in[3]
  T1_2H = development_params.in[4]
  temp_kelvin = temperature.in+273.15
  gas_constant = 8.3144598
  development_rate_timeseries = rho25*(temp_kelvin/298)*exp((DeltaHA/gas_constant)*(1/298-1/temp_kelvin))/(1+exp((DeltaHH/gas_constant)*(1/T1_2H-1/temp_kelvin)))
  return(development_rate_timeseries)
}

#This neglects larval density dependence. Also, we'll use the above (original) function (get_mortality_vector) for adult mortality
calculate_mortality_rate_timeseries = function(mortality_thresholds.in, temperature_max.in){
  T0 = mortality_thresholds.in[1]
  Tinf = mortality_thresholds.in[2]
  survival_rate = array(data = 0.99,dim = c(length(temperature_max.in)))
  temp_mortality_terms = which(temperature_max.in>T0 & temperature_max.in<=Tinf)
  survival_rate[temp_mortality_terms] = 0.99*(1-0.95*(temperature_max.in[temp_mortality_terms]-T0)/(Tinf-T0))
  survival_rate[which(temperature_max.in>Tinf)] = 0.99*0.05
  return(1-survival_rate)
}

calculate_adult_mortality_rate_timeseries = function(mortality_thresholds.in, temperature_max.in){
  T0 = mortality_thresholds.in[1]
  Tinf = mortality_thresholds.in[2]
  survival_rate = array(data = 0.89,dim = c(length(temperature_max.in)))
  temp_mortality_terms = which(temperature_max.in>T0 & temperature_max.in<=Tinf)
  survival_rate[temp_mortality_terms] = 0.89*(1-0.95*(temperature_max.in[temp_mortality_terms]-T0)/(Tinf-T0))
  survival_rate[which(temperature_max.in>Tinf)] = 0.89*0.05
  return(1-survival_rate)
}

differentiate_timeseries = function(series, times){
  #method 1
  # dS1 = dS2 = diff(series)/as.numeric(diff(times))
  # dS1 = append(dS1, dS1[length(dS1)])
  # dS2 = append(dS2, dS2[1], after=0)
  # dS = (dS1+dS2)/2
  #method 2
  t = seq(1, length(times))
  predict(sm.spline(t, series), t, 1)[,1]
}


get_pupal_timeseries = function(adult_surf.in, adult_mortality.in, pupal_dev.in, times.in){
  N = na.spline(adult_surf.in)
  dN = differentiate_timeseries(N, times.in)
  P = (dN+adult_mortality.in*N)/pupal_dev.in
  return(P)
}

get_larval_timeseries = function(pupal_surf.in, pupal_mortality.in, pupal_dev.in, larval_dev.in, times.in){
  dP = differentiate_timeseries(pupal_surf.in, times.in)
  L = (dP+(pupal_mortality.in+pupal_dev.in)*pupal_surf.in)/larval_dev.in
  return(L)
}

#larval inflow refers to egg development * eggs
get_larval_inflow_timeseries = function(larval_surf.in, larval_mortality.in, larval_dev.in, larval_capacity, larval_power, times.in){
  dL = differentiate_timeseries(larval_surf.in, times.in)
  deE = dL + (larval_mortality.in + larval_dev.in + larval_surf.in^larval_power/larval_capacity)*larval_surf.in 
  return(deE)
}

#see methods for this derivation - taking initial larval devlopment from SkeeterBuster.
estimate_initial_eggs = function(egg_mortality0, egg_dev, timestep, adults0, eggs_per_gon_cycle, gonotrophic_cycle_rate0, larval_inflow){
  egg_dev0 = egg_dev[1]
  egg_dev1 = egg_dev[2]
  larval_inflow0 = larval_inflow[1]
  larval_inflow1 = larval_inflow[2]
  eggs0 = (egg_dev0*(0.5*adults0*eggs_per_gon_cycle*gonotrophic_cycle_rate0 - larval_inflow0) - (larval_inflow1-larval_inflow0)/timestep)/(egg_mortality0*egg_dev0-(egg_dev1-egg_dev0)/timestep)
  return(eggs0)
}

egg_model_ode = function(t, x, params) {
  E = x[1]
  N <- params$adult_TS[t]
  gon_rate <- params$gonotrophic_cycle_rate.in[t]
  egg_deaths <- params$egg_mortality.in[t]
  new_larvae <- params$larval_inflow_surf.in[t]
  dE = 0.5*N*params$eggs_per_gon_cycle*gon_rate-egg_deaths*E-new_larvae
  return(list(dE, NULL))
}

egg_model_ode_2 = function(t, x, params) {
  E = x[1]
  N <- params$adult_TS[t]
  gon_rate <- params$gonotrophic_cycle_rate.in[t]
  egg_deaths <- params$egg_mortality.in[t]
  egg_dev <- params$egg_development_rate.in[t]
  dE = 0.5*N*params$eggs_per_gon_cycle*gon_rate-(egg_deaths + egg_dev)*E
  return(list(dE, NULL))
}

larval_model_ode = function(t, x, params) {
  L = x[1]
  P <- params$pupal_surf.in[t]
  E <- params$egg_surf.in[t]
  coeff <- params$larval_coefficient[t]
  egg_dev <- params$egg_development_rate.in[t]
  larval_dev <- params$larval_development_rate.in[t]
  kappa <- params$larval_capacity.in
  alpha <- params$larval_power.in
  dL = egg_dev*E + coeff*L - larval_dev*L^2/P - L^(alpha+1)/kappa
  return(list(dL, NULL))
}


get_egg_timeseries_ode = function(larval_inflow_surf.in, adult_surf.in, egg_mortality.in, gonotrophic_cycle_rate.in, eggs_per_gon_cycle, initial_eggs, times.in){
  times = seq(1,length(times.in),1)
  adult_TS = na.spline(adult_surf.in)
  E = ode(c('eggs' = initial_eggs), times, egg_model_ode, mget(c("eggs_per_gon_cycle", "adult_TS", "gonotrophic_cycle_rate.in", "egg_mortality.in", "larval_inflow_surf.in")))[,2]
  return(E)
}

get_egg_timeseries_ode_2 = function(adult_surf.in, egg_mortality.in, egg_development_rate.in, gonotrophic_cycle_rate.in, eggs_per_gon_cycle, initial_eggs, times.in, ...){
  times = seq(1,length(times.in),1)
  adult_TS = na.spline(adult_surf.in)
  E = ode(c('eggs' = initial_eggs), times, egg_model_ode_2, mget(c("eggs_per_gon_cycle", "adult_TS", "gonotrophic_cycle_rate.in", "egg_mortality.in", "egg_development_rate.in")), ...)[,2]
  return(E)
}

get_larval_timeseries_ode = function(egg_surf.in, pupal_surf.in, egg_development_rate.in, larval_development_rate.in, pupal_development_rate.in, larval_mortality.in, pupal_mortality.in, larval_capacity.in, larval_power.in, initial_larvae, times.in, ...){
  times = seq(1,length(times.in),1)
  dPdt = differentiate_timeseries(pupal_surf.in, times.in)
  larval_coefficient = -larval_mortality.in - larval_development_rate.in + pupal_mortality.in + pupal_development_rate.in + dPdt/pupal_surf.in
  Larv = ode(c('larvae' = initial_larvae), times, larval_model_ode, mget(c("larval_coefficient", "egg_surf.in", "pupal_surf.in", "egg_development_rate.in", "larval_development_rate.in","larval_capacity.in", "larval_power.in")), ...)[,"larvae"]
  return(Larv)
}

get_egg_timeseries = function(larval_surf.in, larval_mortality.in, larval_dev.in, egg_dev.in, larval_capacity, larval_power, times.in){
  dL = differentiate_timeseries(larval_surf.in, times.in)
  E = (dL+(larval_mortality.in + larval_dev.in + larval_surf.in^larval_power/larval_capacity)*larval_surf.in)/egg_dev.in
  return(E)
}

get_female_egg_laying_rate = function(adult_surf.in, egg_surf.in, egg_mortality.in, egg_dev.in, times.in){
  dE = differentiate_timeseries(egg_surf.in, times.in)
  a = (dE+(egg_mortality.in + egg_dev.in)*egg_surf.in)/adult_surf.in
  return(a)
}


get_egg_development_timeseries_alternative_method = function(adult_surf.in, egg_surf.in, eggs_per_gon_cycle, gonotrophic_cycle_rate.in, egg_mortality.in, times.in){
  dE = differentiate_timeseries(egg_surf.in, times.in)
  adult_TS = na.spline(adult_surf.in)
  de = (0.5*adult_TS*eggs_per_gon_cycle*gonotrophic_cycle_rate.in - egg_mortality.in*egg_surf.in - dE)/egg_surf.in
  return(de)
}

full_model_ode = function(t, x, params){
  with(as.list(c(x, params)), {
    devE = params[[1]][t]; devL = params[[2]][t]; devP = params[[3]][t]; devN = params[[4]][t];
    mE = params[[5]][t]; mL = params[[6]][t]; mP = params[[7]][t]; mA = params[[8]][t];
    kL = params[[9]]; alphaL = params[[10]]
    dE = devN*N - (mE + devE)*E
    dL = devE*E - (mL + L^alphaL/kL + devL)*L
    dP = devL*L - (mP + devP)*P
    dN = devP*P - mA*N
    return(list(c(dE, dL, dP, dN)))
  })
}

full_model_ode_2 = function(t, x, params){
  with(as.list(c(x, params)), {
    devE = params[[1]][t]; devL = params[[2]][t]; devP = params[[3]][t]; devN = params[[4]][t];
    mE = params[[5]][t]; mL = params[[6]][t]; mP = params[[7]][t]; mc = params[[8]][t]; mA = params[[9]][t];
    kL = params[[10]]; alphaL = params[[11]]; a = params[[12]]
    dE = devN*a*N/2 - (mE + devE)*E
    dL = devE*E - (mL + L^alphaL/kL + devL + mc)*L
    dP = devL*L - (mP + devP + mc)*P
    dN = devP*P - mA*N
    return(list(c(dE, dL, dP, dN)))
  })
}


Estimate.w.offset.posix = function(XY,tvec,Offset){
  ModelToPlot <- get(sprintf("mod%s",as.character(Offset)))
  EstMat <- StdMat <- matrix(NA, length(XY[,1]),length(tvec))
  XY.sp <- SpatialPoints(XY)
  InorOut <- which(!is.na(over(XY.sp , Boundary)))
  
  TermNames <- c("Multiplier",attr(ModelToPlot$terms,"term.labels"))
  for (tmpTime in 1:length(tvec)){
    DataFrameMat <- matrix(NA,length(XY[,1]),length(TermNames))
    colnames(DataFrameMat) <- TermNames
    
    Location <- which.min(as.numeric(as.POSIXct(Data$Date, tz = "GMT")-tvec[tmpTime])^2)
    
    DataFrameMat[,"Multiplier"] <- Offset
    DataFrameMat[,"YearDay"] <- Data[Location,"YearDay"]
    DataFrameMat[,"Day"] <- Data[Location,"Day"]
    DataFrameMat[,"cWindMax"] <- Data[Location,"cWindMax"]
    DataFrameMat[,"Effort"] <- 0#Data[Location,"Effort"]
    
    DataFrameMat[InorOut,"Xcor"] <- XY[InorOut,1]
    DataFrameMat[InorOut,"Ycor"] <- XY[InorOut,2]
    
    if (length(Include)){
      for (i in Include){
        for (j in 1:CovLen[i]){
          tmpname <- sprintf("%s%i",CovNames[i],j)
          tmpvec <- get(tmpname)
          DataFrameMat[,tmpname] <- tmpvec[Location]
        }
      }
    }
    
    tmp <- predict(ModelToPlot,data.frame(DataFrameMat),se.fit=TRUE)
    EstMat[,tmpTime] <- tmp$fit
    StdMat[,tmpTime] <- tmp$se.fit
  }
  
  out <- vector("list",2)
  out[[1]] <- EstMat
  out[[2]] <- StdMat
  return(out)
}

#=======================================================
# Functions to test each step of the code
#=======================================================
pupal_test_ode = function(t, x, params){
  with(as.list(c(x, params)), {
    devP = params[[1]][t]; 
    P = params[[2]][t]
    mA = params[[3]][t];
    dN = devP*P - mA*N
    return(list(c(dN)))
  })
}

larval_test_ode = function(t, x, params){
  with(as.list(c(x, params)), {
    devL = params[[1]][t]; devP = params[[2]][t];
    L = params[[3]][t]
    mP = params[[4]][t]; mA = params[[5]][t];
    dP = devL*L - mP*P - devP*P
    dN = devP*P - mA*N
    return(list(c(dP, dN)))
  })
}

larval_inflow_test_ode = function(t, x, params){
  with(as.list(c(x, params)), {
    devL = params[[1]][t]; devP = params[[2]][t];
    larval_inflow = params[[3]][t]
    mL = params[[4]][t]; mP = params[[5]][t]; mA = params[[6]][t];
    kL = params[[7]]; alphaL = params[[8]]
    dL = larval_inflow - mL*L - devL*L - L^(1+alphaL)/kL
    dP = devL*L - mP*P - devP*P
    dN = devP*P - mA*N
    return(list(c(dL, dP, dN)))
  })
}

alternative_egg_test_ode = function(t, x, params){
  with(as.list(c(x, params)), {
    devE = params[[1]][t]; devL = params[[2]][t]; devP = params[[3]][t];
    eggs = params[[4]][t]
    mL = params[[5]][t]; mP = params[[6]][t]; mA = params[[7]][t];
    kL = params[[8]]; alphaL = params[[9]]
    dL = eggs*devE - mL*L - devL*L - L^(1+alphaL)/kL
    dP = devL*L - mP*P - devP*P
    dN = devP*P - mA*N
    return(list(c(dL, dP, dN)))
  })
}


