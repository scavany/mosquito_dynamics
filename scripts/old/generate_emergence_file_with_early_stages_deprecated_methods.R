#==============================
#The following tests different methods for calcuating the egg_dev elopment timeseries (no longer used):
#==============================
location = 1
larval_capacity = Inf #set to Inf for no density dependence
adult_timeseries = exp(small_abundance_surf[[1]])[location, ]*proto_pac
pupal_timeseries = get_pupal_timeseries(adult_timeseries, DeathRate, pupal_development_rate, tvec)
pupal_timeseries_nn = pmax(pupal_timeseries, rep(0, length(pupal_timeseries)))
#pupal_timeseries[which(pupal_timeseries<0)]=0
larval_timeseries = get_larval_timeseries(pupal_timeseries, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_timeseries_nn = get_larval_timeseries(pupal_timeseries_nn, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_timeseries_nn = pmax(larval_timeseries_nn, rep(0, length(larval_timeseries_nn)))
#larval_timeseries[which(larval_timeseries<0)]=0

#The next step needs an improved larval mortality (incorporating density dependence)
larval_inflow_timeseries = get_larval_inflow_timeseries(larval_timeseries, larval_mortality, larval_development_rate, larval_capacity, tvec)
larval_inflow_timeseries_nn = get_larval_inflow_timeseries(larval_timeseries_nn, larval_mortality, larval_development_rate, larval_capacity, tvec)
larval_inflow_timeseries_nn = pmax(larval_inflow_timeseries_nn, rep(0, length(larval_inflow_timeseries)))
#larval_inflow_timeseries[which(larval_inflow_timeseries<0)]=0

#estimate the initial eggs, and then the egg timeseries:
eggs_per_gon_cycle = 63 #Otero, 2006
eggs0 = estimate_initial_eggs(egg_mortality[1], egg_development_rate[1:2], as.numeric(diff(tvec))[1], adult_timeseries[1], eggs_per_gon_cycle, gonotrophic_cycle_rate[1], larval_inflow_timeseries[1:2])
eggs0_nn = estimate_initial_eggs(egg_mortality[1], egg_development_rate[1:2], as.numeric(diff(tvec))[1], adult_timeseries[1], eggs_per_gon_cycle, gonotrophic_cycle_rate[1], larval_inflow_timeseries_nn[1:2])
eggs0_nn = max(0, eggs0_nn)
egg_timeseries = get_egg_timeseries(larval_inflow_timeseries, adult_timeseries, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0, tvec)
egg_timeseries_nn = get_egg_timeseries(larval_inflow_timeseries_nn, adult_timeseries, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0_nn, tvec)
egg_timeseries_nn = pmax(egg_timeseries_nn, rep(0, length(egg_timeseries)))

#This is the timeseries we will use as input to the ABM:
egg_development_timeseries_nn1 = larval_inflow_timeseries/egg_timeseries
egg_development_timeseries_nn1 = pmax(egg_development_timeseries_nn1, rep(0, length(egg_development_timeseries_nn1)))
egg_development_timeseries_nn2 = larval_inflow_timeseries_nn/egg_timeseries_nn
egg_development_timeseries_nn2 = pmax(egg_development_timeseries_nn2, rep(0, length(egg_development_timeseries_nn2)))
egg_development_timeseries_nn3 = get_egg_development_timeseries_alternative_method(adult_timeseries, egg_timeseries, eggs_per_gon_cycle, gonotrophic_cycle_rate, egg_mortality, tvec)
egg_development_timeseries_nn3 = pmax(egg_development_timeseries_nn3, rep(0, length(egg_development_timeseries_nn3)))
egg_development_timeseries_nn4 = get_egg_development_timeseries_alternative_method(adult_timeseries, egg_timeseries_nn, eggs_per_gon_cycle, gonotrophic_cycle_rate, egg_mortality, tvec)
egg_development_timeseries_nn4 = pmax(egg_development_timeseries_nn4, rep(0, length(egg_development_timeseries_nn4)))
egg_development_timeseries_5 = larval_inflow_timeseries/egg_timeseries
egg_development_timeseries_6 = get_egg_development_timeseries_alternative_method(adult_timeseries, egg_timeseries, eggs_per_gon_cycle, gonotrophic_cycle_rate, egg_mortality, tvec)

# Methods to compare========================================
# 1. make non-negative at end, use method A
# 2. Make non-negative throughout, use method A
# 3. Make non-negative at end, use method B
# 4. Make non-negative throughout, use method B
# 5. Do not make non-negative, use method A
# 6. Do not make non-negative, use method B
#===========================================================
initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries[1])
initial_conditions_nn = c(E =egg_timeseries_nn[1], L = larval_timeseries_nn[1], P = pupal_timeseries_nn[1], N = adult_timeseries[1])
params_1 = list(egg_development_timeseries_nn1, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity)
params_2 = list(egg_development_timeseries_nn2, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity)
params_3 = list(egg_development_timeseries_nn3, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity)
params_4 = list(egg_development_timeseries_nn4, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity)
params_5 = list(egg_development_timeseries_5, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity)
params_6 = list(egg_development_timeseries_6, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity)
population_timeseries_1 = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_1)
population_timeseries_2 = ode(y = initial_conditions_nn, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_2)
population_timeseries_3 = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_3)
population_timeseries_4 = ode(y = initial_conditions_nn, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_4)
population_timeseries_5 = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_5)
population_timeseries_6 = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params_6)

adult_timeseries_filled = fill_missing_data(adult_timeseries, tvec)
#average_adult_diff = c(mean(abs(adult_timeseries_filled-population_timeseries_1[,"N"]), rm.na=T), mean(abs(adult_timeseries_filled-population_timeseries_2[,"N"]), rm.na=T), mean(abs(adult_timeseries_filled-population_timeseries_3[,"N"]), rm.na=T),mean(abs(adult_timeseries_filled-population_timeseries_4[,"N"]), rm.na=T), mean(abs(adult_timeseries_filled-population_timeseries_5[,"N"]), rm.na=T), mean(abs(adult_timeseries_filled-population_timeseries_6[,"N"]), rm.na=T))
average_adult_diff = c(mean(abs(adult_timeseries_filled-population_timeseries_1[,"N"]), rm.na=T), mean(abs(adult_timeseries_filled-population_timeseries_2[,"N"]), rm.na=T), mean(abs(adult_timeseries_filled-population_timeseries_5[,"N"]), rm.na=T), mean(abs(adult_timeseries_filled-population_timeseries_6[,"N"]), rm.na=T))

#Plot these
n=length(adult_timeseries)
plotvals = seq(1, n, 24)
plot(adult_timeseries[plotvals])
par(mfrow=c(3,2))
plot(population_timeseries_1[plotvals,"N"])
plot(population_timeseries_2[plotvals,"N"])
#plot(population_timeseries_3[plotvals,"N"])
#plot(population_timeseries_4[plotvals,"N"])
plot(population_timeseries_5[plotvals,"N"])
plot(population_timeseries_6[plotvals,"N"])

#Plot timeseries 5 with data
par(mfrow=c(1,2))
plot(adult_timeseries[plotvals])
plot(population_timeseries_5[plotvals,"N"])


#===================================================================
# Estimate the timeseries in chunks - deprecated_method
#===================================================================
# Note that there are 4018 time points, which is divisible by 7, so do 7 chunks
location = 1
larval_capacity = Inf #set to Inf for no density dependence
eggs_per_gon_cycle = 63 #Otero, 2006
number_chunks = 100
chunk_size = length(tvec)/number_chunks
egg_timeseries = c()
egg_development_timeseries =c()
adult_timeseries = exp(small_abundance_surf[[1]])[location,]*proto_pac
pupal_timeseries = get_pupal_timeseries(adult_timeseries, DeathRate, pupal_development_rate, tvec)
larval_timeseries = get_larval_timeseries(pupal_timeseries, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_inflow_timeseries = get_larval_inflow_timeseries(larval_timeseries, larval_mortality, larval_development_rate, larval_capacity, tvec)
for (i in 1:number_chunks){
  start_pt = floor((i-1)*chunk_size+1 + 0.5)
  end_pt = floor(i*chunk_size + 0.5)
  #For subsequent steps, better to use previous chunk to estimate eggs0?
  eggs0 = estimate_initial_eggs(egg_mortality[start_pt], egg_development_rate[start_pt:(start_pt+1)], as.numeric(diff(tvec))[start_pt], adult_timeseries[start_pt], eggs_per_gon_cycle, gonotrophic_cycle_rate[start_pt], larval_inflow_timeseries[start_pt:(start_pt+1)])
  egg_timeseries = c(egg_timeseries, get_egg_timeseries(larval_inflow_timeseries[start_pt:end_pt], adult_timeseries[start_pt:end_pt], egg_mortality[start_pt:end_pt], gonotrophic_cycle_rate[start_pt:end_pt], eggs_per_gon_cycle, eggs0, tvec[start_pt:end_pt]))
  egg_development_timeseries = c(egg_development_timeseries, larval_inflow_timeseries[start_pt:end_pt]/egg_timeseries[start_pt:end_pt])
} 

initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries[1])
params = list(egg_development_timeseries, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity) 
population_timeseries = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params)

# That didn't seem to work, but it looks like the pattern is correct, but the magnitude grows? So multiply by a factor?

#===========================================================================
# Try and find a correction function, which is a function of time.
# Try f1 = linear
# and f2 = quadratic, both =1 when t=0
# and f3 = exponential
#===========================================================================
location = 1
larval_capacity = Inf #set to Inf for no density dependence
eggs_per_gon_cycle = 63 #Otero, 2006
adult_timeseries = exp(small_abundance_surf[[1]])[location,]*proto_pac
pupal_timeseries = get_pupal_timeseries(adult_timeseries, DeathRate, pupal_development_rate, tvec)
larval_timeseries = get_larval_timeseries(pupal_timeseries, pupal_mortality, pupal_development_rate, larval_development_rate, tvec)
larval_inflow_timeseries = get_larval_inflow_timeseries(larval_timeseries, larval_mortality, larval_development_rate, larval_capacity, tvec)
eggs0 = estimate_initial_eggs(egg_mortality[1], egg_development_rate[1:2], as.numeric(diff(tvec))[1], adult_timeseries[1], eggs_per_gon_cycle, gonotrophic_cycle_rate[1], larval_inflow_timeseries[1:2])
egg_timeseries = get_egg_timeseries(larval_inflow_timeseries, adult_timeseries, egg_mortality, gonotrophic_cycle_rate, eggs_per_gon_cycle, eggs0, tvec)
egg_development_timeseries = larval_inflow_timeseries/egg_timeseries
initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries[1])
params = list(egg_development_timeseries, larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity) 
population_timeseries_corrected = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params)
initial_conditions = c(E =egg_timeseries[1], L = larval_timeseries[1], P = pupal_timeseries[1], N = adult_timeseries[1])

plot(population_timeseries[plotvals,"N"]/adult_timeseries[plotvals])
#Relationship looks exponential - try a exponential growth in egg_development

times_minus_1 = 1:length(tvec) - 1

k_series = log(fill_missing_data(adult_timeseries, tvec)[2:length(tvec)]/population_timeseries[2:length(tvec),"N"])/times_minus_1[2:length(tvec)]
k_guess = mean(k_series)
for (i in seq(10)){
  k_guess = k_guess/2
  params = list(egg_development_timeseries*exp(k_guess*times_minus_1), larval_development_rate, pupal_development_rate, gonotrophic_cycle_rate, eggs_per_gon_cycle, egg_mortality, larval_mortality, pupal_mortality, DeathRate, larval_capacity) 
  population_timeseries_corrected = ode(y = initial_conditions, times = seq(1,length(tvec),1), func = full_model_ode, parms = params)
  plot(population_timeseries_corrected[,"N"])
}


adult_timeseries_filled = fill_missing_data(adult_timeseries, tvec)
mean_errors = mean(abs(adult_timeseries_filled-population_timeseries_5[,"N"])/adult_timeseries_filled)

