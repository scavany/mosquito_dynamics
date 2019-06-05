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

#original method
pop.model.ode <- function (t, x, params) {
  N <- x[1];
  cor_factor <- params[1]
  mu <- DeathRate[t+1]; epsilon <- Emergence_TS[t] ;
  dN <-  epsilon * cor_factor - mu * N
  return(list(dN))
}

# #either of new methods
# pop.model.ode <- function (t, x, params) {
#   N <- x[1]; 
#   cor_factor <- params[1]
#   mu <- DeathRate[t]; 
#   epsilon <- Emergence.data$Emergence_Timeseries[t] ; 
#   dN <-  epsilon*cor_factor - mu * N
#   return(list(dN))
# }