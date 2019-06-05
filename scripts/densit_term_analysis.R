load("surfaces.RData")
times=seq(length(tvec))

#alpha and beta taken from Legros et al. 2009, using 63 eggs taken from Otero et al.
# these are daily rates, so need to account for the gonotrophic-cycle rate``
beta = 0.302
alpha = (0.229)*(63*mean(gonotrophic_cycle_rate)^beta)

burn_period = 180
dvec = seq(as.Date("2000-01-01")-burn_period,as.Date("2010-12-31")+burn_period,by="1 day")
tvec = seq(as.Date("2000-01-01")-burn_period,as.Date("2010-12-31")+burn_period,by="1 day")

alpha_list = 10^seq(-10,2, 0.1)
submin_vector = vector(length=length(alpha_list))
min_vector = vector(length=length(alpha_list))
max_vector = vector(length=length(alpha_list))
median_vector = vector(length=length(alpha_list))
for (i in seq(length(alpha_list))){
  submin_vector[i] = exp(-alpha_list[i]*0.01^beta)
  min_vector[i] = exp(-alpha_list[i]*0.1^beta)
  max_vector[i] = exp(-alpha_list[i]*2.6^beta)
  median_vector[i] = exp(-alpha_list[i]*1.0^beta)
}

plot(alpha_list, submin_vector/median_vector, type='l', log='x')

beta_list = 10^seq(-10,0.2, 0.1)
submin_vector = vector(length=length(beta_list))
min_vector = vector(length=length(beta_list))
max_vector = vector(length=length(beta_list))
median_vector = vector(length=length(beta_list))
for (i in seq(length(beta_list))){
  submin_vector[i] = exp(-alpha*0.01^beta_list[i])
  min_vector[i] = exp(-alpha*0.1^beta_list[i])
  max_vector[i] = exp(-alpha*2.6^beta_list[i])
  median_vector[i] = exp(-alpha*1.0^beta_list[i])
}

plot(beta_list, submin_vector/median_vector, type='l', log='x')

