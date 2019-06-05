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
  dN = differentiate_timeseries(adult_surf.in, times.in)
  P = (dN+adult_mortality.in*adult_surf.in)/pupal_dev.in
  return(P)
}

get_larval_timeseries = function(pupal_surf.in, pupal_mortality.in, pupal_dev.in, larval_dev.in, times.in){
  dP = differentiate_timeseries(pupal_surf.in, times.in)
  L = (dP+(pupal_mortality.in+pupal_dev.in)*pupal_surf.in)/larval_dev.in
  return(L)
}

egg_model_ode = function(t, x, params) {
  E = x[1]
  N <- params$adult_TS[t]
  adult_dev <- params$adult_development_rate.in[t]
  egg_deaths <- params$egg_mortality.in[t]
  egg_dev <- params$egg_development.in[t]
  dE = N*adult_dev-egg_deaths*E-egg_dev*E
  return(list(dE, NULL))
}

get_egg_timeseries = function(adult_surf.in, egg_mortality.in, adult_development_rate.in, egg_development.in, initial_eggs, times.in){
  times = seq(1,length(times.in),1)
  adult_TS = adult_surf.in
  E = ode(c('eggs' = initial_eggs), times, egg_model_ode, mget(c("adult_TS", "egg_mortality.in", "adult_development_rate.in", "egg_development.in")))[,2]
  return(E)
}

get_larval_capacity_timeseries = function(larval_surf.in, egg_surf.in, larval_mortality.in, larval_development_rate.in, egg_development_rate.in, larval_power.in, times.in) {
  dL = differentiate_timeseries(larval_surf.in, times.in)
  kl = (as.vector(larval_surf.in)^(1+larval_power))/(egg_development_rate.in*egg_surf.in - larval_mortality.in*larval_surf.in - larval_development_rate.in*larval_surf.in - dL)
  return(kl)
}

full_model_ode = function(t, x, params){
  with(as.list(c(x, params)), {
    devE = params[[8]][t]; devL = params[[7]][t]; devP = params[[6]][t]; devN = params[[5]][t];
    mE = params[[4]][t]; mL = params[[3]][t]; mP = params[[2]][t]; mA = params[[1]][t];
    K = params[[9]][t]; alpha = params[[10]]
    dE = devN*N - (mE + devE)*E
    dL = devE*E - (mL + L^alpha/K + devL)*L
    dP = devL*L - (mP + devP)*P
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
    N=x[1]
    devP = params[[3]][t]; 
    P = params[[1]][t]
    mA = params[[2]][t];
    dN = devP*P - mA*N
    return(list(c(dN)))
  })
}

larval_test_ode = function(t, x, params){
  with(as.list(c(x, params)), {
    P = x[1];
    devL = params[[4]][t]; devP = params[[3]][t];
    L = params[[1]][t]
    mP = params[[2]][t]; 
    dP = devL*L - mP*P - devP*P
    return(list(c(dP)))
  })
}

#=======================================================
# Integration functions
#=======================================================
evalN = function(x, A){
  answer = array(0,length(x))
  answer[which(x<1)] = A[1]
  answer[which(x>=1)] = A[x[which(x>=1)]]
  return(answer)
}

emergence_integrand_mean = function(n, t, adults, lambdas.in){
  return(evalN(t-n, adults)*dhypoexp(n, lambdas.in))
}

density_integrand_varying = function(n, t, egg_dev, larval_dev, pupal_dev){
  lambdas = matrix(c(evalN(t-n, egg_dev), evalN(t-n, larval_dev), evalN(t-n, pupal_dev)), ncol=3)
  density_vector = vector(length=length(t))
  for (i in seq(length(t))){
    density_vector[i] = dhypoexp(n, lambdas[i,])
  }
  return(density_vector)
}

emergence_integrand_varying = function(n, t, adults, egg_dev, larval_dev, pupal_dev){
  return(evalN(t-n, adults)*density_integrand_varying(n, t, egg_dev, larval_dev, pupal_dev))
}

integrate_trapezoid = function(integrand, lower, upper, nstrips,  ...){
  Dx = (upper-lower)/nstrips
  sum_of_strips = 0
  for (i in seq(nstrips)){
    strip = integrand((i-1)*Dx,...) + integrand(i*Dx,...)
    sum_of_strips = sum_of_strips + strip
  }
  return(sum_of_strips*Dx/2)
}

density_hypo_exp3 = function(x, l1, l2, l3){
  return(l1*l2*l3*((l3-l1)*exp(-x*l2) + (l1-l2)*exp(-x*l3) + (l2-l3)*exp(-x*l1))/((l1-l3)*(l2-l1)*(l3-l2)))
}

Lambda = function(t, rate_vector){
  answer = matrix(data=0,nrow=length(t), ncol=length(t))
  for (i in seq(2, length(t))){
    for (j in seq(i-1)){
      answer[i,j] = (sum(rate_vector[j:(i-1)]) + sum(rate_vector[(j+1):i]))/2
    }
  }
  return(answer)
}

dnonhomexp = function(t, rate_vector, Lambda_){
  answer = matrix(data=0,nrow=length(t), ncol=length(t))
  for (i in seq(1, length(t))){
    for (j in seq(i)){
      answer[i,j] = rate_vector[i]*exp(-Lambda_[i,j])
    }
  }
  return(answer)
}

convolution = function(f_1, f_2, nmax){
  if (dim(f_1)[1]!=dim(f_2)[1] | dim(f_1)[2]!=dim(f_2)[2]){
    return(print("matrices different shape"))
  }
  maxi = dim(f_1)[1] 
  maxj = dim(f_1)[2] 
  answer = matrix(0, maxi, maxj)
  #Do j loop first
  for (j in seq(maxj)){
    print(j)
    for (i in seq(j, min(nmax+j,maxi))){
      f_left = c(vector(length=maxi+j-2),f_1[j:maxi,j])
      f_right = c(vector(length=i+j-2), rev(f_2[j:maxi,j]), vector(length=maxi-i))
      f_prod = f_left*f_right
      answer[i,j] = (sum(f_prod[1:(length(f_prod)-1)]) + sum(f_prod[2:length(f_prod)]))/2
    }
  }
  return(answer)
}

